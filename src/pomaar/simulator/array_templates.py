"""
Polarization-aware unit-cell templates for the array builder.

An array element is built from a unit-cell template (an HFSS design) whose port sheets carry
their native polarization in the name (`PortSheet_H`, `PortSheet1_V`). Placing a template at a
yaw rotates its native E-field: the built polarization of a port is R(yaw) applied to its native
axis, which must lie along x (H) or y (V) of the global frame. This module plans the elements
(template choice, active and load ports, names, signs) and checks them before anything is drawn.

It is pure Python/NumPy and imports nothing from AEDT, so it is unit-tested on the host.
"""

from dataclasses import dataclass, field
import math
import re

POLARIZATIONS = ("H", "V")
# Global polarization basis: H along x, V along y (Ludwig-3 X/Y in the array PhaseCentreCS)
POLARIZATION_AXES = {"H": (1.0, 0.0), "V": (0.0, 1.0)}
# Matched wave ports on the unused ports of a dual-port cell end in this suffix. SBR+ port
# parsing rejects them (no trailing H/V), so they are never Tx or Rx.
LOAD_SUFFIX = "_Load"
# Separator between a template object's name and its template key in the target design
TEMPLATE_SEPARATOR = "__"
POLARIZATION_TOLERANCE_DEG = 1.0


class LayoutError(ValueError):
    """Raised when the layout and the unit-cell templates are inconsistent."""


class PolarizationMismatchError(LayoutError):
    """Raised when a placed template does not build the polarization its element declares."""


@dataclass(frozen=True)
class PortSpec:
    """A port sheet of a unit-cell template."""

    sheet: str
    polarization: str | None  # native polarization 'H'/'V' (None if untagged and undeclared)
    phase_centre: tuple  # (x, y, z) in mm, template coordinates


@dataclass
class UnitCellTemplate:
    """
    A unit-cell design used as an element template. The pure planning code uses `design`, `key`,
    `ports` and `offset`; the builder fills in the AEDT handles and the object classification.
    """

    design: str
    ports: tuple = ()
    offset: tuple = (0.0, 0.0, 0.0)  # element reference point (PhaseCentreCS), template mm
    key: str = ""
    app: object = None
    active: list = field(default_factory=list)
    port_sheets: list = field(default_factory=list)
    dummies: dict = field(default_factory=dict)  # (operation, target) -> [object names]
    layers: dict = field(default_factory=dict)  # name -> {material, z_min, z_max}

    def __post_init__(self):
        if not self.key:
            self.key = template_key(self.design)

    @property
    def native_polarizations(self):
        """The template's tagged polarizations, e.g. ('V',) or ('H', 'V')."""
        return tuple(p for p in POLARIZATIONS if any(s.polarization == p for s in self.ports))

    @property
    def objects(self):
        """All objects replicated per element: active copper, port sheets and dummy solids."""
        dummies = [name for names in self.dummies.values() for name in names]
        return list(self.active) + list(self.port_sheets) + dummies

    def tagged(self, name):
        """The object's name in the target design while it serves as a copied template."""
        return f"{name}{TEMPLATE_SEPARATOR}{self.key}"


@dataclass
class PlannedPort:
    """A port of a placed element."""

    spec: PortSpec
    name: str  # wave-port name in the array design
    sheet: str  # replicated sheet name in the array design
    polarization: str | None  # built polarization (None if not along x or y)
    sign: int  # +1 if E points along +x (H) / +y (V), -1 if opposite, 0 if unknown
    axis: tuple  # built E-field unit vector (x, y) in global coordinates
    active: bool  # False for a matched load on an unused port
    cs_name: str | None  # far-field phase-reference CS (None for loads and dropped ports)
    phase_centre: tuple  # (x, y, z) in mm, array coordinates


@dataclass
class ElementPlan:
    """An element of the layout with its template, placement and ports."""

    label: str
    raw_label: str
    role: str | None
    position: tuple  # (x, y, z) in mm: where the template's offset point lands
    position_expression: tuple | None  # parametric (x, y[, z]) or None
    yaw: float
    polarizations: tuple  # requested polarizations, e.g. ('V',) or ('H', 'V'); () = unchecked
    template: UnitCellTemplate
    ports: list = field(default_factory=list)


def template_key(design_name):
    """An AEDT-safe object-name tag for a template design."""
    return re.sub(r"\W", "_", design_name)


def port_sheet_polarization(sheet_name):
    """The polarization tag of a port sheet ('PortSheet_V' -> 'V', 'PortSheet1_h' -> 'H')."""
    match = re.search(r"_([HhVv])$", sheet_name)
    return match.group(1).upper() if match else None


def parse_polarizations(value):
    """Parses an element's polarization ('v', 'hv', ['h', 'v'], None) into e.g. ('H', 'V')."""
    if value is None:
        return ()
    letters = "".join(value) if isinstance(value, (list, tuple)) else str(value)
    letters = letters.replace(",", "").replace(" ", "").replace("_", "").upper()
    if any(letter not in POLARIZATIONS for letter in letters):
        raise LayoutError(f"Unknown polarization '{value}'; use h, v or hv.")
    return tuple(p for p in POLARIZATIONS if p in letters)


def unit_cell_key(polarizations):
    """The `--unit-cell` / `metadata.unit_cells` key of a polarization set: 'H', 'V' or 'HV'."""
    return "".join(polarizations)


def parse_unit_cell_args(values):
    """
    Parses `--unit-cell` values into {key: design}: 'h=PatchH' and 'v=PatchV' serve single-pol
    elements, 'hv=DualPatch' serves dual-polarized elements (polarization: hv).
    """
    mapping = {}
    for value in values or []:
        key, sep, design = str(value).partition("=")
        if not sep or not design.strip():
            raise LayoutError(
                f"--unit-cell expects '<pol>=<Design>' (e.g. v=PatchV), got '{value}'."
            )
        polarizations = parse_polarizations(key.strip())
        if not polarizations:
            raise LayoutError(
                f"--unit-cell needs a polarization (h, v or hv) before '=', got '{value}'."
            )
        _assign_unit_cell(mapping, unit_cell_key(polarizations), design.strip(), "--unit-cell")
    return mapping


def merge_unit_cells(cli_mapping, metadata_mapping):
    """
    Merges the `--unit-cell` mapping with a layout's `metadata.unit_cells` ({h: ..., v: ...,
    hv: ...}). The same key mapped to two different designs is an error.
    """
    merged = dict(cli_mapping or {})
    for key, design in (metadata_mapping or {}).items():
        polarizations = parse_polarizations(key)
        if polarizations:
            _assign_unit_cell(
                merged, unit_cell_key(polarizations), str(design), "metadata.unit_cells"
            )
    return merged


def _assign_unit_cell(mapping, key, design, source):
    if mapping.get(key, design) != design:
        raise LayoutError(
            f"{source} maps polarization {key.lower()} to '{design}', but it is already mapped "
            f"to '{mapping[key]}'."
        )
    mapping[key] = design


def rotate_xy(vector, yaw_deg):
    """Rotates (x, y) by yaw about the z axis."""
    c, s = math.cos(math.radians(yaw_deg)), math.sin(math.radians(yaw_deg))
    return (c * vector[0] - s * vector[1], s * vector[0] + c * vector[1])


def place_xy(point, offset, yaw_deg, position):
    """
    Lateral placement rule of the builder, p -> R(yaw) (p - offset) + position: the template's
    offset point (its PhaseCentreCS) lands on the layout position.
    """
    dx, dy = rotate_xy((point[0] - offset[0], point[1] - offset[1]), yaw_deg)
    return (dx + position[0], dy + position[1])


def place(point, offset, z_ref, yaw_deg, position):
    """
    Full placement: lateral by `place_xy`, vertical by the common stack shift z -> z - z_ref +
    position_z, so that every template shares one board stack.
    """
    x, y = place_xy(point, offset, yaw_deg, position)
    z_pos = position[2] if len(position) > 2 else 0.0
    return (x, y, point[2] - z_ref + z_pos)


def built_polarization(native, yaw_deg, tolerance_deg=POLARIZATION_TOLERANCE_DEG):
    """
    The polarization a port of native polarization `native` radiates when its template is
    placed at `yaw_deg`: returns (pol, sign, axis) with pol 'H'/'V' or None when the rotated
    axis is not within `tolerance_deg` of x or y, sign +1/-1 along +x/+y, and the axis itself.
    """
    if native not in POLARIZATION_AXES:
        return None, 0, (0.0, 0.0)
    axis = rotate_xy(POLARIZATION_AXES[native], yaw_deg)
    limit = math.cos(math.radians(tolerance_deg))
    for pol, reference in POLARIZATION_AXES.items():
        projection = axis[0] * reference[0] + axis[1] * reference[1]
        if abs(projection) >= limit:
            return pol, 1 if projection > 0 else -1, axis
    return None, 0, axis


def resolve_template(polarizations, explicit, unit_cells, templates, primary):
    """
    Picks an element's template: its explicit `unit_cell`, else the unit cell mapped to its
    polarization set ('h', 'v' or 'hv'; a dual-polarized element also takes a design mapped to
    both 'h' and 'v'), else the primary template (polarization by rotation, checked later).
    """
    if explicit:
        if explicit not in templates:
            raise LayoutError(f"Unknown unit cell '{explicit}' (loaded: {', '.join(templates)}).")
        return templates[explicit]
    key = unit_cell_key(polarizations)
    if key in unit_cells:
        return templates[unit_cells[key]]
    designs = {unit_cells.get(pol) for pol in polarizations}
    if len(polarizations) == 2 and len(designs) == 1 and None not in designs:
        return templates[designs.pop()]
    return templates[primary]


def strip_polarization_suffix(label, pol):
    """'Rx_1_V' -> 'Rx_1' for pol 'V'; other labels unchanged."""
    suffix = f"_{pol}"
    return label[: -len(suffix)] if pol and label.upper().endswith(suffix) else label


def element_label(raw_label, polarizations):
    """A single-polarized element's label carries its polarization ('Rx_1' -> 'Rx_1_V')."""
    if len(polarizations) == 1:
        pol = polarizations[0]
        return raw_label if raw_label.upper().endswith(f"_{pol}") else f"{raw_label}_{pol}"
    return raw_label


def port_name(label, pol, index=None, load=False):
    """
    Wave-port name: 'Port_<label>[_F<n>]_<P>[_Load]', with the element label stripped of its own
    polarization suffix, so that a single-pol 'Rx_1_V' gives 'Port_Rx_1_V' and the two ports of
    a dual-pol 'Rx_1' give 'Port_Rx_1_H' / 'Port_Rx_1_V' (SBR+ pairs them by the shared prefix).
    """
    core = strip_polarization_suffix(label, pol)
    feed = f"_F{index}" if index else ""
    pol_part = f"_{pol}" if pol else ""
    return f"Port_{core}{feed}{pol_part}{LOAD_SUFFIX if load else ''}"


def cs_name_for_port(port):
    """Element coordinate system of a port: 'Port_Rx_1_V' -> 'CS_Rx_1_V'."""
    return f"CS_{port[len('Port_') :]}" if port.startswith("Port_") else f"CS_{port}"


def parse_port(port_name):
    """
    Classifies a port name as (role, polarization, element key), e.g. 'TxH' -> ('Tx', 'H', 'Tx')
    and 'Port_Tx_1_H' -> ('Tx', 'H', 'Port_Tx_1_'); None if it is neither (e.g. a '_Load' port).
    The element key pairs the H and V ports of one dual-polarized element. Excitation names of a
    linked antenna ('<instance>_Port_Tx_1_H') parse the same way.
    """
    roles = re.findall(r"Tx|Rx", port_name)
    pol = re.search(r"([HV])$", port_name)
    if not roles or not pol:
        return None
    # The last Tx/Rx counts, so that a prefix such as a linked antenna's instance name
    # ('TxRxDemo1_Port_Rx_1_H') cannot change the role
    return roles[-1], pol.group(1), port_name[:-1]


def is_load_port(port_name):
    """True for the matched load on an unused port of a dual-port cell."""
    return port_name.split(":")[0].endswith(LOAD_SUFFIX)


def plan_element(element, templates, primary, unit_cells, drop_unused_ports=False, z_ref=0.0):
    """
    Plans one layout element: template, label, and every port with its built polarization,
    sign, name and phase centre. Raises `PolarizationMismatchError` when the placed template
    cannot build the declared polarization(s).
    """
    raw_label = str(element.get("label", element.get("name", "Element")))
    polarizations = parse_polarizations(element.get("polarization", element.get("pol")))
    position = [float(v) for v in element.get("position", element.get("pos", [0.0, 0.0, 0.0]))]
    position = tuple(position + [0.0] * (3 - len(position)))
    expression = element.get("position_expression", element.get("pos_expr"))
    yaw = float(element.get("yaw", element.get("rotation_yaw", element.get("rotation", 0.0))))
    template = resolve_template(
        polarizations, element.get("unit_cell"), unit_cells, templates, primary
    )
    label = element_label(raw_label, polarizations)

    role = element.get("role")
    if role is not None and str(role) not in ("Tx", "Rx"):
        raise LayoutError(f"Element '{label}': role must be Tx or Rx, not '{role}'.")

    if not template.ports:
        raise LayoutError(f"Unit cell '{template.design}' has no port sheet.")

    built = [(spec, *built_polarization(spec.polarization, yaw)) for spec in template.ports]
    if polarizations:
        missing = [p for p in polarizations if not any(b[1] == p for b in built)]
        if missing:
            raise PolarizationMismatchError(
                _mismatch_message(label, polarizations, template, yaw, built, unit_cells)
            )

    # Repeated polarizations within one element get a feed index
    counts = {}
    for spec, pol, _, _ in built:
        key = pol or spec.polarization
        counts[key] = counts.get(key, 0) + 1
    # Port names start from the label without the element's own polarization suffix
    core = strip_polarization_suffix(label, polarizations[0]) if len(polarizations) == 1 else label
    seen = {}
    ports = []
    for spec, pol, sign, axis in built:
        active = not polarizations or pol in polarizations
        name_pol = pol or spec.polarization
        seen[name_pol] = seen.get(name_pol, 0) + 1
        index = seen[name_pol] if counts[name_pol] > 1 else None
        name = port_name(core, name_pol, index, load=not active)
        if not active and drop_unused_ports:
            name = None
        ports.append(
            PlannedPort(
                spec=spec,
                name=name,
                sheet=f"{spec.sheet}_{label}",
                polarization=pol,
                sign=sign,
                axis=axis,
                active=active,
                cs_name=cs_name_for_port(name) if active else None,
                phase_centre=place(spec.phase_centre, template.offset, z_ref, yaw, position),
            )
        )

    for port in ports:
        if not port.active or not port.polarization:
            continue
        parsed = parse_port(port.name)
        if parsed is not None and role is None:
            role = parsed[0]
        if parsed is None:
            raise LayoutError(
                f"Port '{port.name}' of element '{label}' has no Tx/Rx in its name, so SBR+ "
                "cannot use it; put Tx or Rx in the element label."
            )
        if role is not None and parsed[0] != role:
            raise LayoutError(
                f"Port '{port.name}' reads as {parsed[0]}, but element '{label}' has role {role}."
            )

    return ElementPlan(
        label=label,
        raw_label=raw_label,
        role=role,
        position=position,
        position_expression=tuple(str(e) for e in expression) if expression else None,
        yaw=yaw,
        polarizations=polarizations,
        template=template,
        ports=ports,
    )


def _mismatch_message(label, polarizations, template, yaw, built, unit_cells):
    wanted = "".join(polarizations).lower()
    natives = "/".join(template.native_polarizations) or "untagged"
    results = ", ".join(f"{spec.sheet} -> {pol or 'neither H nor V'}" for spec, pol, _, _ in built)
    hints = []
    if any(spec.polarization is None for spec, _, _, _ in built):
        hints.append(
            "tag its port sheets (PortSheet_H / PortSheet_V) or declare their native "
            "polarization (--native-polarization for the source design, --unit-cell for others)"
        )
    elif len(polarizations) == 2 and len(template.native_polarizations) < 2:
        hints.append(
            "a dual-polarized element needs a dual-port unit cell (--unit-cell hv=<Design>)"
        )
    else:
        absent = [
            p
            for p in polarizations
            if p not in unit_cells and p not in template.native_polarizations
        ]
        turn = "turn the element by 90 deg"
        if absent:
            article = "an" if absent[0] == "H" else "a"
            hints.append(
                f"add {article} {''.join(absent)}-polarized unit cell "
                f"(--unit-cell {''.join(absent).lower()}=<Design>) or {turn}"
            )
        else:
            hints.append(turn)
    return (
        f"Element '{label}' declares polarization {wanted}, but unit cell '{template.design}' "
        f"({natives}-native) at yaw {yaw:g} deg builds {results}; " + " ".join(hints) + "."
    )


def plan_elements(elements, templates, primary, unit_cells, drop_unused_ports=False, z_ref=0.0):
    """
    Plans all layout elements and checks the array as a whole: unique labels and port names.
    Returns (plans, warnings).
    """
    plans = [
        plan_element(element, templates, primary, unit_cells, drop_unused_ports, z_ref)
        for element in elements
    ]
    warnings = []

    labels = [plan.label for plan in plans]
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise LayoutError(f"Duplicate element labels: {', '.join(duplicates)}.")

    names = [port.name for plan in plans for port in plan.ports if port.name]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise LayoutError(f"Duplicate port names: {', '.join(duplicates)}.")

    for i, first in enumerate(plans):
        for second in plans[i + 1 :]:
            if math.dist(first.position[:2], second.position[:2]) < 1e-6:
                warnings.append(
                    f"Elements '{first.label}' and '{second.label}' share the position "
                    f"{list(first.position[:2])} mm, so their geometries overlap; use one "
                    "dual-polarized element (polarization: hv) with a dual-port unit cell."
                )
    return plans, warnings


def map_duplicates(pasted_names, tagged_to_base):
    """
    Maps the names AEDT gives duplicated template objects back to their template object, e.g.
    'PortSheet1__PatchV_1' -> 'PortSheet1'. A duplicate is named after its source plus an
    optional '_' and a number; exactly one source must match, and every source exactly once.
    """
    mapping = {}
    for pasted in pasted_names:
        candidates = [
            tagged
            for tagged in tagged_to_base
            if pasted == tagged or re.fullmatch(re.escape(tagged) + r"_?\d+", pasted)
        ]
        if len(candidates) != 1:
            raise RuntimeError(
                f"Cannot map duplicated object '{pasted}' to one template object "
                f"(candidates: {candidates or 'none'})."
            )
        mapping[pasted] = tagged_to_base[candidates[0]]
    bases = list(mapping.values())
    repeated = sorted({base for base in bases if bases.count(base) > 1})
    missing = sorted(set(tagged_to_base.values()) - set(bases))
    if repeated or missing:
        raise RuntimeError(
            f"Duplicated template objects do not match the templates (repeated: {repeated}, "
            f"missing: {missing})."
        )
    return mapping


def offset_expression(expression, delta_mm):
    """A parametric coordinate plus a numeric offset in mm (omitted below 0.1 um)."""
    if abs(delta_mm) < 1e-4:
        return str(expression)
    return f"({expression}) + ({delta_mm:.4f}mm)"


def port_table(plans):
    """One row per port (active, load and dropped) for printing and the JSON sidecar."""
    rows = []
    for plan in plans:
        for port in plan.ports:
            rows.append(
                {
                    "port": port.name,
                    "element": plan.label,
                    "role": plan.role,
                    "polarization": port.polarization,
                    "sign": port.sign,
                    "e_axis": [round(v, 6) for v in port.axis],
                    "active": port.active,
                    "template": plan.template.design,
                    "template_sheet": port.spec.sheet,
                    "native_polarization": port.spec.polarization,
                    "yaw_deg": plan.yaw,
                    "phase_centre_mm": [round(v, 6) for v in port.phase_centre],
                }
            )
    return rows


def pair_signs(rows):
    """
    The pair sign pi = s_H * s_V of each role, or None where a role lacks H or V ports or mixes
    signs within one polarization. The measured matrix is diag(1, pi_Rx) S diag(1, pi_Tx) up to
    a common sign: pi_Rx = -1 negates the receive-V row (S_VH, S_VV), pi_Tx = -1 the transmit-V
    column (S_HV, S_VV).
    """
    signs = {}
    for row in rows:
        if row["active"] and row["role"] and row["polarization"] and row["sign"]:
            signs.setdefault(row["role"], {}).setdefault(row["polarization"], set()).add(
                row["sign"]
            )
    result = {}
    for role in ("Tx", "Rx"):
        by_pol = signs.get(role, {})
        if all(len(by_pol.get(pol, ())) == 1 for pol in POLARIZATIONS):
            result[role] = next(iter(by_pol["H"])) * next(iter(by_pol["V"]))
        else:
            result[role] = None
    return result


def sign_warnings(rows):
    """Warnings about field signs that a sphere calibration cannot reveal."""
    warnings = []
    by_role_pol = {}
    for row in rows:
        if row["active"] and row["role"] and row["polarization"] and row["sign"]:
            key = (row["role"], row["polarization"])
            by_role_pol.setdefault(key, {}).setdefault(row["sign"], []).append(row["port"])
    for (role, pol), groups in sorted(by_role_pol.items()):
        if len(groups) > 1:
            warnings.append(
                f"{role} {pol} ports radiate with opposite signs (+: {', '.join(groups[1])}; "
                f"-: {', '.join(groups[-1])}); their channels differ by 180 deg."
            )
    signs = pair_signs(rows)
    flipped = {"Rx": "receive-V row (S_VH, S_VV)", "Tx": "transmit-V column (S_HV, S_VV)"}
    for role in ("Tx", "Rx"):
        if signs.get(role) == -1:
            warnings.append(
                f"{role} H and V ports have opposite signs (pi_{role} = -1): the measured "
                f"{flipped[role]} is negated relative to S_HH. A sphere calibration cannot "
                "reveal this; change the yaw of the H or V elements by 180 deg."
            )
    return warnings


def sidecar_dict(design_name, plans, z_ref):
    """Contents of the '<design>_ports.json' sidecar."""
    rows = port_table(plans)
    return {
        "design": design_name,
        "convention": (
            "Global frame: H along +x, V along +y. sign = +1 if a port's E-field at boresight "
            "points along +x (H) or +y (V). Measured S = diag(s_RxH, s_RxV) S "
            "diag(s_TxH, s_TxV)."
        ),
        "z_reference_mm": z_ref,
        "pair_signs": pair_signs(rows),
        "ports": rows,
    }


# ----------------------------------------------------------------------------------------------
# Channels of the linked array (shared by sbr_simulator, sbr_results and polarimetry_processor)
# ----------------------------------------------------------------------------------------------

QUAD_KEYS = ("RxH", "RxV", "TxH", "TxV")


def s_parameter_ports(expression):
    """'S(<receive>,<transmit>)' -> (receive port, transmit port); None for other expressions."""
    match = re.fullmatch(r"\s*S\(\s*([^,()]+?)\s*,\s*([^,()]+?)\s*\)\s*", str(expression))
    return (match[1], match[2]) if match else None


def _element_label(key, role):
    """Short label of an element key: 'Linked1_Port_Rx_1_' -> 'Rx_1', 'Horn1_Rx' -> 'Rx'."""
    start = key.rfind(role)
    return key[start:].strip("_") if start >= 0 else key.strip("_")


def polarimetric_quads(ports):
    """
    Groups Tx/Rx port (or excitation) names into polarimetric quads, {name: {'RxH', 'RxV',
    'TxH', 'TxV': port}}, one per pair of a dual-polarized Rx and a dual-polarized Tx. Per role,
    the dual-polarized elements are the H/V port pairs that share an element key (`parse_port`);
    a role without any uses its single H and V ports when it has exactly one of each (e.g. the
    polarimetric cluster, whose H and V elements sit apart). Quad names read '<Rx>/<Tx>'.
    Load ports and names without Tx/Rx and H/V are ignored.
    """
    keyed = {"Tx": {}, "Rx": {}}
    single = {"Tx": {"H": [], "V": []}, "Rx": {"H": [], "V": []}}
    for port in dict.fromkeys(ports):
        parsed = None if is_load_port(port) else parse_port(port)
        if parsed is None:
            continue
        role, pol, key = parsed
        keyed[role].setdefault(key, {})[pol] = port
        single[role][pol].append(port)

    elements = {}
    for role in ("Rx", "Tx"):
        dual = {
            _element_label(key, role): pols
            for key, pols in keyed[role].items()
            if set(pols) == set(POLARIZATIONS)
        }
        if not dual and all(len(single[role][pol]) == 1 for pol in POLARIZATIONS):
            dual = {role: {pol: single[role][pol][0] for pol in POLARIZATIONS}}
        elements[role] = dual

    return {
        f"{rx_name}/{tx_name}": {
            "RxH": rx["H"],
            "RxV": rx["V"],
            "TxH": tx["H"],
            "TxV": tx["V"],
        }
        for rx_name, rx in elements["Rx"].items()
        for tx_name, tx in elements["Tx"].items()
    }


def select_quad(ports, quad=None):
    """
    Resolves a quad: a dict {'RxH', 'RxV', 'TxH', 'TxV': port} is returned as it is, a name is
    looked up in `polarimetric_quads(ports)`, and None takes the only quad there is. Returns
    (name, quad).
    """
    if isinstance(quad, dict):
        missing = [key for key in QUAD_KEYS if key not in quad]
        if missing:
            raise ValueError(f"Quad needs the ports {', '.join(QUAD_KEYS)}; missing {missing}.")
        return "custom", {key: quad[key] for key in QUAD_KEYS}
    quads = polarimetric_quads(ports)
    if quad is not None:
        if quad not in quads:
            raise ValueError(f"Unknown quad '{quad}'; available: {', '.join(quads) or 'none'}.")
        return quad, quads[quad]
    if len(quads) != 1:
        hint = (
            f"choose one of {', '.join(quads)}"
            if quads
            else "no dual-polarized Tx/Rx pair was found; pass the four ports explicitly as "
            "{'RxH': ..., 'RxV': ..., 'TxH': ..., 'TxV': ...}"
        )
        raise ValueError(f"Cannot pick a polarimetric quad automatically: {hint}.")
    return next(iter(quads.items()))


def attach_port_table(table, excitations):
    """
    Matches a builder port table (`sidecar_dict`) to the excitation names of a linked antenna
    ('<instance>_<port>', where the instance is the array design's name plus a number, e.g.
    'PatchCluster1_Port_Tx_H'). Returns the table's rows that match, each with its 'excitation',
    or None unless every Tx/Rx excitation (loads excepted) is found.
    """
    pattern = re.compile(re.escape(table["design"]) + r"\d*_(.+)")
    by_port = {row["port"]: row for row in table["ports"] if row.get("port")}
    rows = []
    for excitation in dict.fromkeys(excitations):
        match = pattern.fullmatch(excitation)
        row = by_port.get(match[1]) if match else None
        if row is None:
            if is_load_port(excitation) or parse_port(excitation) is None:
                continue
            return None
        rows.append({**row, "excitation": excitation})
    return rows or None
