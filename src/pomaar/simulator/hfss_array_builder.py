#!/usr/bin/env python3
"""
Module 1: HFSS Full-Wave Array Synthesis.
Programmatically builds a planar MIMO antenna array layout on a single PCB substrate in HFSS
using a template-based, boolean-driven assembly workflow.
"""

import json
import math
import os
import re
import sys

from ansys.aedt.core import Desktop, Hfss
from ansys.aedt.core.generic.constants import Axis, Gravity
from ansys.aedt.core.generic.data_handlers import _dict2arg
from ansys.aedt.core.generic.general_methods import grpc_active_sessions
from ansys.aedt.core.internal.load_aedt_file import load_entire_aedt_file
import numpy as np

from pomaar.simulator.array_templates import (
    LayoutError,
    PortSpec,
    UnitCellTemplate,
    array_design_name,
    is_load_port,
    layout_settings,
    map_duplicates,
    offset_expression,
    place,
    place_xy,
    plan_elements,
    port_sheet_polarization,
    port_table,
    sidecar_dict,
    sign_warnings,
)
from pomaar.simulator.frequency_sweeps import (
    LINKABLE_SWEEP_TYPE,
    adaptive_coverage_warning,
    adaptive_frequencies_ghz,
    design_setups,
    native_sweep_arguments,
    setup_sweeps,
    sweep_bounds_ghz,
    sweep_point_count,
    sweep_ranges,
)
from pomaar.simulator.frequency_sweeps import parse_frequency_ghz as _parse_frequency_ghz

SPEED_OF_LIGHT_MM_S = 2.99792458e11
DEFAULT_CENTRE_FREQUENCY_GHZ = 79.0
# The array's analysis setup (the unit cell's meshing) and the builder's discrete sweep in it
ARRAY_SETUP_NAME = "Setup"
ARRAY_SWEEP_NAME = "Sweep"

# Builder-owned far-field spheres. They are (re)created on every results setup so that a
# manual GUI edit (e.g. a changed Boresight) can never silently alter the reported quantities.
FF_SPHERE_3D = "InfiniteSphere"
FF_SPHERE_PHI0 = "Phi0 cut"
FF_SPHERE_PHI90 = "Phi90 cut"

# Builder-owned far-field reports. They depend on the sphere definitions above and are
# therefore recreated together with them.
FF_REPORT_3D = "Realized gain 3D"
FF_REPORT_PHI0 = "Realized gain Phi0"
FF_REPORT_PHI90 = "Realized gain Phi90"
FF_REPORT_XPD_PHI0 = "XPD Phi0"
FF_REPORT_XPD_PHI90 = "XPD Phi90"

# Pattern quantities. Realized gain is used for absolute patterns. XPD is a co/cross ratio and is
# therefore identical for Directivity, Gain and Realized gain (they differ only by a common
# power normalization); directivity is used so that XPD does not depend on port mismatch settings.
# With 'Linear' polarization (slant 0) and an unrotated CS, CoPolar/CrossPolar are Ludwig-3 X/Y
# referenced to the sphere CS x-axis: XPD of an x-polarized (H) port is Co - Cross, and of a
# y-polarized (V) port Cross - Co.
XPD_X_POL_EXPRESSION = "dB(DirCoPolar) - dB(DirCrossPolar)"
XPD_Y_POL_EXPRESSION = "dB(DirCrossPolar) - dB(DirCoPolar)"


class BuildAbortedError(Exception):
    """Raised when the user declines to continue at an interactive prompt."""


class NativeAedtRefusedError(RuntimeError):
    """Raised when connecting from a Linux host outside the AEDT container."""


class NoLiveSessionError(RuntimeError):
    """Raised when no AEDT GUI session is running to attach to."""


def connect_live_session():
    """
    Attaches to the running AEDT GUI session, the pipeline's only mode of operation, and returns
    its Desktop. Nothing here ever starts an AEDT.

    On a Linux host PyAEDT cannot see a containerized AEDT's gRPC socket (separate network
    namespace) and would silently launch a native AEDT instead, which crashes there: run the
    scripts inside the container (`ansys-vnc exec ...`), or set POMAAR_ALLOW_NATIVE_AEDT=1 on a
    native installation.
    """
    if (
        sys.platform.startswith("linux")
        and not os.path.exists("/run/.containerenv")
        and os.environ.get("POMAAR_ALLOW_NATIVE_AEDT") != "1"
    ):
        raise NativeAedtRefusedError(
            "Refusing to connect from outside the AEDT container: run this via "
            "`ansys-vnc exec`, or set POMAAR_ALLOW_NATIVE_AEDT=1 for a native AEDT install."
        )
    ports = grpc_active_sessions(non_graphical=False)
    if not ports:
        raise NoLiveSessionError(
            "No AEDT GUI session is running; start one (`ansys-vnc aedt`) and open the project."
        )
    if len(ports) > 1:
        print(f"Several AEDT GUI sessions are running (gRPC ports {ports}); using {ports[0]}.")
    print(f"Attaching to the AEDT GUI session on gRPC port {ports[0]}...")
    return Desktop(port=ports[0], new_desktop=False, non_graphical=True, close_on_exit=False)


def _confirm(question, preset=None, default=True):
    """
    Resolves a yes/no decision: a preset (from CLI flags) wins, then an interactive prompt,
    then the default when stdin is not a terminal.
    """
    if preset is not None:
        return bool(preset)
    if not sys.stdin.isatty():
        return default
    hint = "y" if default else "n"
    answer = input(f"{question} (y/n) [default: {hint}]: ").strip().lower()
    if not answer:
        return default
    return answer not in ["n", "no"]


def _setup_frequency_ghz(app):
    """Returns the representative (centre) frequency of the first analysis setup in GHz, or None."""
    setup_names = list(app.setup_names)
    if not setup_names:
        return None
    setup = app.get_setup(getattr(app, "active_setup", None) or setup_names[0])
    # The single adaptive frequency, else the mean of the multi-frequency or broadband ones
    adaptive = adaptive_frequencies_ghz(setup.props)
    if adaptive:
        return sum(adaptive) / len(adaptive)

    for sweep_name in setup.get_sweep_names():
        sweep_props = getattr(setup.get_sweep(sweep_name), "props", {})
        start = _parse_frequency_ghz(sweep_props.get("RangeStart"))
        stop = _parse_frequency_ghz(sweep_props.get("RangeEnd"))
        if start and stop:
            return (start + stop) / 2.0
    return None


def saved_setups(app):
    """
    {setup: properties, each with its 'Sweeps'} of the design of `app`, read from the project
    file after saving the project, so that the file holds the current definitions. PyAEDT's
    setup and sweep objects are not used: PyAEDT 1.1 sets every loaded sweep's 'Type' to
    'Interpolating' whatever it is, and reads a project file that may predate GUI edits.
    """
    app.save_project()
    return design_setups(load_entire_aedt_file(app.project_file), app.design_name)


def _length_mm(app, expression):
    """
    Converts a length expression to mm. Plain literals are parsed directly so that PyAEDT does not
    create its temporary 'pyaedt_evaluator' variable in the design.
    """
    match = re.fullmatch(
        r"\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*(mm|um|nm|cm|m|metres?|meters?|mil|in)?\s*",
        str(expression),
    )
    if match:
        unit = match.group(2)
        # AEDT writes metres as 'meter'; British and plural spellings are accepted as well
        if unit and unit.startswith(("metre", "meter")):
            unit = "m"
        scale = {None: 1.0, "mm": 1.0, "um": 1e-3, "nm": 1e-6, "cm": 10.0, "m": 1e3,
                 "mil": 0.0254, "in": 25.4}
        return float(match.group(1)) * scale[unit]
    return float(app.evaluate_expression(expression)) * 1000.0


def _spelling_variants(name):
    """
    The name itself plus its British/American counterparts (Centre/Center, Metre/Meter), British
    first, so that designs, coordinate systems and variables named in either spelling are found.
    """
    pairs = (("Centre", "Center"), ("centre", "center"), ("Metre", "Meter"), ("metre", "meter"))
    british, american = name, name
    for british_part, american_part in pairs:
        british = british.replace(american_part, british_part)
        american = american.replace(british_part, american_part)
    return list(dict.fromkeys([name, british, american]))


def _find_by_spelling(names, name):
    """Returns the first of `name`'s spelling variants present in `names`, or None."""
    return next((variant for variant in _spelling_variants(name) if variant in names), None)


def _delete_coordinate_system(modeler, cs_name):
    """Deletes a coordinate system by name if it exists."""
    for cs in list(modeler.coordinate_systems):
        if cs.name == cs_name:
            try:
                cs.delete()
            except Exception:
                pass


def _delete_far_field_sphere(app, sphere_name):
    """Deletes a far-field infinite sphere setup by name if it exists."""
    for setup in list(app.field_setups):
        if setup.name == sphere_name:
            setup.delete()


def _z_extent(obj):
    """Returns (z_min, z_max) of an object's bounding box."""
    bbox = obj.bounding_box
    z_coords = [float(bbox[2]), float(bbox[5])]
    return min(z_coords), max(z_coords)


def _far_field_columns(solution_data, expressions):
    """
    Unpacks a PyAEDT far-field SolutionData into flat arrays: one per intrinsic sweep variable
    (e.g. 'Theta', 'Phi', 'Freq') and one complex array per expression.
    """
    real, imag = solution_data.full_matrix_real_imag
    intrinsic_names = list(solution_data.intrinsics.keys())
    first = np.asarray(real[expressions[0]])
    offset = first.shape[1] - 1 - len(intrinsic_names)  # leading design-variation columns
    columns = {name: first[:, offset + i] for i, name in enumerate(intrinsic_names)}
    for expr in expressions:
        columns[expr] = np.asarray(real[expr])[:, -1] + 1j * np.asarray(imag[expr])[:, -1]
    return columns


class MimoHfssBuilder:
    """
    Automates the creation of a planar MIMO array on a single PCB substrate
    in HFSS full-wave from an isolated antenna element design by using naming
    conventions and boolean operations on replicated dummy solids.
    """

    def __init__(self, project_path, layout, target_design_name=None):
        """
        `layout` is the path of a layout YAML/JSON file or its loaded dict. The layout is bound to
        its unit cells (see `array_templates.layout_settings`): its `primary_unit_cell` is the
        source design, which also provides the analysis setup, and the array design is named
        `<primary_unit_cell><LayoutFileStem>` (`<primary_unit_cell>MimoArray` for a dict).
        """
        layout_path = None
        if not isinstance(layout, dict):
            layout_path = os.path.abspath(layout)
            layout = load_layout_file(layout_path)
        self.layout = layout
        self.layout_path = layout_path
        settings = layout_settings(layout)
        self.project_path = os.path.abspath(project_path)
        self.source_design_name = settings.primary_unit_cell
        if target_design_name:
            self.target_design_name = target_design_name
        elif layout_path:
            self.target_design_name = array_design_name(settings.primary_unit_cell, layout_path)
        else:
            self.target_design_name = f"{settings.primary_unit_cell}MimoArray"
        # Taken from the primary unit cell's setup once it is open, like the frequency sweep
        self.centre_frequency_ghz = None
        self.phase_reference_mode = settings.phase_reference_mode
        self.unit_cells = settings.unit_cells
        self.drop_unused_ports = settings.drop_unused_ports
        self.port_to_cs_mapping = {}
        self.element_phase_centres = []

        self.desktop_session = None
        self.source_design_app = None
        self.target_design_app = None
        self.is_new_design = True

    def connect_desktop(self):
        """Attaches to the running AEDT GUI session (`connect_live_session`)."""
        self.desktop_session = connect_live_session()

    def _open_design(self, design_name, project=None):
        """Opens (and activates) a design of the project as its own Hfss app."""
        return Hfss(
            project=project or self.project_path,
            design=design_name,
            new_desktop=False,
            close_on_exit=False,
        )

    def _classify_template(self, app, design_name):
        """
        STEP 1 for one unit cell: classifies its objects by name convention into global layers,
        port sheets (tagged with their polarization, e.g. 'PortSheet_V'), active copper and
        boolean dummy solids.
        """
        template = UnitCellTemplate(design=design_name, app=app)
        print(f"\nScanning and classifying objects in unit cell '{design_name}':")
        for name in list(app.modeler.object_names):
            obj = app.modeler.get_object_from_name(name)
            if not obj:
                continue

            # Check for boolean dummy solids first (e.g. Subtract_L2_Ground)
            if name.startswith("Subtract_") or name.startswith("Unite_"):
                tokens = name.split("_")
                key = (tokens[0], "_".join(tokens[1:]))
                template.dummies.setdefault(key, []).append(name)
                print(f"  [Dummy {key[0]}]   {name} -> target: {key[1]}")
                continue

            # Skip vacuum domain solids, radiation boundaries, or airboxes present in source design
            if (
                obj.material_name.lower() == "vacuum"
                or "RadiatingSurface" in name
                or "Airbox" in name
                or "RadiationBox" in name
            ):
                print(f"  [Skipped Boundary] {name} ({obj.material_name})")
                continue

            is_substrate = name.startswith(("L12_", "L23_", "L34_", "L45_"))
            is_ground = name.startswith(("L2_Ground", "Ground_Plane", "GND"))
            if is_substrate or is_ground:
                z_min, z_max = _z_extent(obj)
                template.layers[name] = {
                    "material": obj.material_name,
                    "z_min": z_min,
                    "z_max": z_max,
                }
                tag = "[Global Layer]    " if is_substrate else "[Global Ground]   "
                print(
                    f"  {tag} {name} ({obj.material_name}, Z=[{z_min:.2f}, {z_max:.2f}] mm)"
                )

            elif "PortSheet" in name or name.startswith("Port_"):
                template.port_sheets.append(name)
                pol = port_sheet_polarization(name)
                print(f"  [Port Sheet]       {name} (polarization {pol or 'untagged'})")

            else:
                template.active.append(name)
                print(f"  [Active Element]   {name}")

        # Port specs get their phase centres in `_template_phase_centres`
        template.ports = tuple(
            PortSpec(sheet, port_sheet_polarization(sheet), (0.0, 0.0, 0.0))
            for sheet in template.port_sheets
        )
        return template

    def _template_phase_centres(self, template, use_existing_cs=None):
        """
        Sets the template's reference point (its `PhaseCentreCS`, reused or fitted) and the phase
        centre of each port. A dual-port cell may hold `PhaseCentreCS_H` / `PhaseCentreCS_V` for
        its ports; it is never fitted, since that needs per-port sources in the unit cell.
        """
        app = template.app
        cs_map = {cs.name: cs for cs in app.modeler.coordinate_systems}
        # Also finds a unit cell's 'PhaseCenterCS'
        main_cs = _find_by_spelling(cs_map, "PhaseCentreCS")
        port_cs = {
            pol: _find_by_spelling(cs_map, f"PhaseCentreCS_{pol}")
            for pol in template.native_polarizations
        }
        port_cs = {pol: name for pol, name in port_cs.items() if name}
        multi_port = len(template.ports) > 1

        def origin(cs_name):
            props = cs_map[cs_name].props
            try:
                return [
                    round(_length_mm(app, props[key]), 3)
                    for key in ("OriginX", "OriginY", "OriginZ")
                ]
            except Exception as e:
                print(f"  Warning: Could not evaluate {cs_name} origin ({e}); using zero offset.")
                return [0.0, 0.0, 0.0]

        offset = [0.0, 0.0, 0.0]
        run_fit = False
        if main_cs:
            offset = origin(main_cs)
            print(f"\n[INFO] {main_cs} found in '{template.design}'. Origin offset: {offset} mm.")
            if not multi_port and not _confirm(
                f"Do you want to use the existing {main_cs}?", preset=use_existing_cs
            ):
                run_fit = True
        elif port_cs:
            centres = [origin(name) for name in port_cs.values()]
            offset = [round(sum(c[i] for c in centres) / len(centres), 3) for i in range(3)]
            print(
                f"\n[INFO] '{template.design}' has {', '.join(port_cs.values())}; reference point "
                f"(their mean): {offset} mm."
            )
        else:
            print("\n" + "=" * 80)
            print(f"[WARNING] PhaseCentreCS was NOT found in unit cell '{template.design}'!")
            print("=" * 80 + "\n")
            if use_existing_cs:
                print("[INFO] Proceeding without PhaseCentreCS offset [0.0, 0.0, 0.0] mm.")
            elif multi_port:
                print(
                    "[INFO] A dual-port unit cell is not fitted (it would need per-port "
                    "sources); add PhaseCentreCS or PhaseCentreCS_H/_V to it."
                )
                if not _confirm("Proceed without PhaseCentreCS offset?"):
                    raise BuildAbortedError("Aborted by user: no phase centre available.")
            elif _confirm("Do you want to fit the phase centre from the unit-cell far field?"):
                run_fit = True
            elif _confirm("Proceed without PhaseCentreCS offset?"):
                print("[INFO] Proceeding with zero offset [0.0, 0.0, 0.0] mm.")
            else:
                raise BuildAbortedError("Aborted by user: no phase centre available.")

        if run_fit:
            offset = self.fit_unit_cell_phase_centre(app=app)

        template.offset = tuple(offset)
        template.ports = tuple(
            PortSpec(
                spec.sheet,
                spec.polarization,
                tuple(origin(port_cs[spec.polarization]))
                if spec.polarization in port_cs
                else template.offset,
            )
            for spec in template.ports
        )

    def _load_templates(self, unit_cells, explicit_designs=(), use_existing_cs=None):
        """
        Opens and classifies the primary unit cell, every design of `unit_cells` ({pol: design})
        and the designs named by elements (`unit_cell:`), and determines their phase centres.
        Returns {design: template}, primary first.
        """
        designs = list(
            dict.fromkeys([self.source_design_name, *unit_cells.values(), *explicit_designs])
        )
        missing = [d for d in designs if d not in self.source_design_app.design_list]
        if missing:
            raise LayoutError(f"Unit-cell design(s) not found in the project: {', '.join(missing)}")

        templates = {}
        for design in designs:
            app = (
                self.source_design_app
                if design == self.source_design_name
                else self._open_design(design, self.source_design_app.project_name)
            )
            template = self._classify_template(app, design)
            self._template_phase_centres(template, use_existing_cs)
            templates[design] = template
        return templates

    @staticmethod
    def _check_shared_stack(templates):
        """All templates must share the board stack (layer names, materials and z extents)."""
        primary, *others = templates.values()
        reference = {
            name: (info["material"], round(info["z_min"], 4), round(info["z_max"], 4))
            for name, info in primary.layers.items()
        }
        for template in others:
            stack = {
                name: (info["material"], round(info["z_min"], 4), round(info["z_max"], 4))
                for name, info in template.layers.items()
            }
            if stack != reference:
                raise LayoutError(
                    f"Unit cell '{template.design}' has a different board stack than "
                    f"'{primary.design}':\n  {stack}\n  vs {reference}"
                )

    @staticmethod
    def _merged_design_variables(templates):
        """
        Design variables of all templates, primary first. A name with different expressions in
        two templates is an error, since the array design has only one value for it.
        """
        merged, sources, conflicts = {}, {}, []
        for template in templates.values():
            for name, variable in template.app.variable_manager.design_variables.items():
                expression = str(getattr(variable, "expression", variable)).replace(" ", "")
                if name in merged and sources[name][1] != expression:
                    conflicts.append(
                        f"{name}: '{sources[name][1]}' in '{sources[name][0]}', "
                        f"'{expression}' in '{template.design}'"
                    )
                    continue
                merged.setdefault(name, variable)
                sources.setdefault(name, (template.design, expression))
        if conflicts:
            raise LayoutError(
                "Unit cells define the same design variable differently; rename variables per "
                "template:\n  " + "\n  ".join(conflicts)
            )
        return merged

    def synthesize_array_in_hfss(self, use_existing_cs=None, overwrite=None, simulate=None):
        """
        Synthesizes the builder's layout in HFSS using a naming-convention boolean assembly:
        1. Classifies the objects of every unit-cell template and plans the elements (template,
           ports, polarization check) before the target design is touched.
        2. Creates or clears the target design.
        3. Sizes and draws the global continuous substrate and ground layers.
        4. Copies each template into the target design exactly ONCE and replicates it to its
           elements with the local duplicate command.
        5. Applies boolean operations using the dummy solids (e.g. Subtract_L2_Ground).
        6. Assigns wave ports (lumped fallback), including matched loads on unused feeds.
        7. Creates the radiation airbox, verifies the placement, reports the port polarizations,
           sets up the simulation, the far-field spheres and the reports, and optionally solves.

        `use_existing_cs`, `overwrite` and `simulate` preset the prompts for reusing the unit
        cells' PhaseCentreCS, replacing an existing array design and solving (`run_solve`).

        Raises
        ------
        BuildAbortedError
            If the user declines to continue at one of the interactive prompts.
        LayoutError
            If the layout and the unit cells are inconsistent (raised before any change).
        """
        layout_metadata = self.layout.get("metadata") or {}
        elements_list = self.layout["elements"]

        if not self.desktop_session:
            self.connect_desktop()

        for key in layout_settings(self.layout).ignored_keys:
            print(f"Note: layout metadata '{key}' is no longer used and is ignored.")

        print(f"Loading project: {self.project_path}")
        self.source_design_app = self._open_design(self.source_design_name)

        # Centre frequency of the primary unit cell's setup (the array copies that setup): it
        # sets the Airbox clearance, the far-field report frequency and a phase-centre fit
        self.centre_frequency_ghz = (
            _setup_frequency_ghz(self.source_design_app) or DEFAULT_CENTRE_FREQUENCY_GHZ
        )
        operating_frequency_ghz = self.centre_frequency_ghz
        print(
            f"Centre frequency: {operating_frequency_ghz:.3f} GHz "
            f"(setup of '{self.source_design_name}')"
        )

        # --- STEP 1: Classify the unit-cell templates and plan the elements ---
        # Everything that can reject the layout runs here, before Step 2 clears a design.
        unit_cells = self.unit_cells
        explicit_designs = [str(e["unit_cell"]) for e in elements_list if e.get("unit_cell")]
        templates = self._load_templates(unit_cells, explicit_designs, use_existing_cs)
        primary = templates[self.source_design_name]
        if not list(primary.app.setup_names):
            raise LayoutError(
                f"Unit cell '{self.source_design_name}' has no analysis setup; the array copies "
                "its meshing and frequency sweep from it."
            )
        self._check_shared_stack(templates)
        design_variables = self._merged_design_variables(templates)
        # All templates share the board stack, which is shifted by the primary unit cell's
        # phase-centre height: z_ref ends up at the layout z of every element.
        z_ref = primary.offset[2]

        plans, layout_warnings = plan_elements(
            elements_list,
            templates,
            self.source_design_name,
            unit_cells,
            drop_unused_ports=self.drop_unused_ports,
            z_ref=z_ref,
        )
        used_templates = list({id(p.template): p.template for p in plans}.values())
        print(f"\nPlanned {len(plans)} elements:")
        for plan in plans:
            ports = ", ".join(
                f"{p.name or p.spec.sheet + ' (dropped)'}" for p in plan.ports
            )
            print(f"  {plan.label}: '{plan.template.design}' at yaw {plan.yaw:g} deg -> {ports}")
        for warning in layout_warnings:
            print(f"  [WARNING] {warning}")

        project_name = self.source_design_app.project_name
        design_list = self.source_design_app.design_list

        # Save project to flush any recently deleted designs from memory cache
        try:
            self.source_design_app.save_project()
        except Exception:
            pass

        # --- STEP 2: Reuse/Create Array Design ---
        # To avoid the buggy delete-and-recreate race condition in AEDT,
        # we check if the design exists. If so, we reuse it and clear its modeler.
        if self.target_design_name in design_list:
            if not _confirm(
                f"Design '{self.target_design_name}' already exists. "
                "Do you want to clear it and overwrite?",
                preset=overwrite,
                default=False,
            ):
                raise BuildAbortedError(
                    f"Kept the existing array design '{self.target_design_name}' "
                    "(--overwrite replaces it)."
                )

            print(f"Design '{self.target_design_name}' already exists. Reusing and clearing it...")
            self.target_design_app = self._open_design(self.target_design_name, project_name)
            self.is_new_design = False

            # Clear existing 3D bodies
            target_modeler = self.target_design_app.modeler
            all_objs = list(target_modeler.object_names)
            if all_objs:
                print(f"  Clearing {len(all_objs)} existing objects...")
                target_modeler.delete(all_objs)

            # Clear coordinate systems
            target_modeler.set_working_coordinate_system("Global")
            for cs in list(target_modeler.coordinate_systems):
                try:
                    cs.delete()
                except Exception:
                    pass

            # Clear any leftover Optimetrics setups to prevent validation errors with deleted setups
            try:
                opt_names = list(self.target_design_app.ooptimetrics.GetSetupNames())
                if opt_names:
                    print(f"  Clearing {len(opt_names)} existing Optimetrics setups...")
                    self.target_design_app.ooptimetrics.DeleteSetups(opt_names)
            except Exception:
                pass

            # Setup and sweep purging will be handled at the end of synthesis in configure_simulation_setup
        else:
            print(f"Creating new HFSS array design '{self.target_design_name}'...")
            self.is_new_design = True
            try:
                self.target_design_app = Hfss(
                    project=project_name,
                    design=self.target_design_name,
                    solution_type="Modal",
                    new_desktop=False,
                    close_on_exit=False,
                )
            except Exception as e:
                print(f"  Failed to create design ({e}). Attempting to connect to cached design...")
                self.target_design_app = self._open_design(self.target_design_name, project_name)

        # Reconnect the template apps so that they point to the unit-cell designs again
        for template in reversed(list(templates.values())):
            template.app = self._open_design(template.design, project_name)
        self.source_design_app = primary.app

        # Sync design variables from the unit cells to the target design
        # This prevents copied objects in the target design from retaining stale/cached design variable dimensions.
        try:
            print("Syncing design variables from the unit cells to the target design...")
            for var_name, var_expr in design_variables.items():
                self.target_design_app[var_name] = var_expr
        except Exception as e:
            print(f"  Warning: Failed to sync design variables ({e})")

        target_modeler = self.target_design_app.modeler

        # Measure the unit-cell reach relative to each template's PhaseCentreCS (offset), from
        # the active elements and port sheets so that the measurement reflects the true element
        # reach (L_feed) and is never corrupted by oversized substrate boundaries.
        uc_x_extents = []
        uc_y_extents = []
        template_corners = {}
        for template in used_templates:
            modeler = template.app.modeler
            offset = template.offset
            for obj_name in list(template.layers) + template.active + template.port_sheets:
                try:
                    bb = modeler.get_object_from_name(obj_name).bounding_box
                    uc_x_extents.extend(
                        [abs(float(bb[0]) - offset[0]), abs(float(bb[3]) - offset[0])]
                    )
                    uc_y_extents.extend(
                        [abs(float(bb[1]) - offset[1]), abs(float(bb[4]) - offset[1])]
                    )
                except Exception:
                    pass
            template_corners[id(template)] = self._template_corners(template)

        unit_cell_extent_x = max(uc_x_extents) if uc_x_extents else 5.0
        unit_cell_extent_y = max(uc_y_extents) if uc_y_extents else 5.0
        unit_cell_extent_max = max(unit_cell_extent_x, unit_cell_extent_y)

        # Calculate clearance using lambda_0 / 4 (quarter wavelength) rule of thumb
        wavelength_mm = SPEED_OF_LIGHT_MM_S / (operating_frequency_ghz * 1e9)
        airbox_clearance_mm = 0.25 * wavelength_mm

        # Inject geometric constants into target design
        self.target_design_app["unitCellExtentX"] = f"{unit_cell_extent_x:.4f}mm"
        self.target_design_app["unitCellExtentY"] = f"{unit_cell_extent_y:.4f}mm"
        self.target_design_app["unitCellExtent"] = f"{unit_cell_extent_max:.4f}mm"
        self.target_design_app["airboxClearance"] = f"{airbox_clearance_mm:.4f}mm"

        # Inject layout variables from layout metadata
        layout_vars = layout_metadata.get("variables", {})
        for var_name, var_expr in layout_vars.items():
            self.target_design_app[var_name] = str(var_expr)

        # --- STEP 3: Size & Draw Global Contiguous Board ---
        # Sort layers by prefix (L1, L12, L2...) to ensure they are created in stackup order
        global_layers = primary.layers
        sorted_layers = sorted(global_layers.keys(), key=lambda name: name.split("_")[0])
        ground_layer_names = []

        # The board covers every element's template footprint (its layer outline in x, the
        # reach of its feedlines and ports in y), placed like the element itself
        all_x_glob, all_y_glob = [], []
        for plan in plans:
            for corner in template_corners[id(plan.template)]:
                x, y = place_xy(corner, plan.template.offset, plan.yaw, plan.position)
                all_x_glob.append(x)
                all_y_glob.append(y)

        # Determine arrayBoardWidth and arrayBoardLength formulas or values
        # We use arrayBoardWidth / arrayBoardLength for the synthesized board layers so that
        # any internal unit-cell variables (like boardWidth / boardLength / sub_L) retain their
        # single-element values without double-expanding the feedlines. The layout's formulas
        # (with its own margin variable, e.g. pcbMargin) win; the tight footprint is a fallback.
        board_config = layout_metadata.get("board", {})
        width_formula = board_config.get("width_formula")
        length_formula = board_config.get("length_formula")

        if all_x_glob:
            board_width_calc = max(all_x_glob) - min(all_x_glob)
            board_height_calc = max(all_y_glob) - min(all_y_glob)
        else:
            board_width_calc = board_height_calc = 10.0

        if width_formula:
            self.target_design_app["arrayBoardWidth"] = str(width_formula)
        else:
            self.target_design_app["arrayBoardWidth"] = f"{board_width_calc:.4f}mm"

        if length_formula:
            self.target_design_app["arrayBoardLength"] = str(length_formula)
        else:
            self.target_design_app["arrayBoardLength"] = f"{board_height_calc:.4f}mm"

        print("\nConstructing global contiguous board layers:")
        for layer_name in sorted_layers:
            layer_info = global_layers[layer_name]
            z_min = layer_info["z_min"]
            z_max = layer_info["z_max"]
            thickness = z_max - z_min
            material = layer_info["material"]

            if layer_name.endswith("_Substrate"):
                print(
                    f"  Creating substrate '{layer_name}' (material={material}, thickness={thickness:.2f} mm)"
                )
                sub_board = target_modeler.create_box(
                    origin=[
                        "-arrayBoardWidth / 2",
                        "-arrayBoardLength / 2",
                        f"{z_min - z_ref:.4f}mm",
                    ],
                    sizes=["arrayBoardWidth", "arrayBoardLength", f"{thickness:.4f}mm"],
                    name=layer_name,
                    material=material,
                )
                sub_board.transparency = 0.5

            elif layer_name.endswith("_Ground"):
                ground_layer_names.append(layer_name)
                if thickness == 0.0:
                    print(
                        f"  Creating ground plane sheet '{layer_name}' at Z={z_min - z_ref:.2f} mm"
                    )
                    target_modeler.create_rectangle(
                        orientation="XY",
                        origin=[
                            "-arrayBoardWidth / 2",
                            "-arrayBoardLength / 2",
                            f"{z_min - z_ref:.4f}mm",
                        ],
                        sizes=["arrayBoardWidth", "arrayBoardLength"],
                        name=layer_name,
                    )
                else:
                    print(
                        f"  Creating ground plane block '{layer_name}' (thickness={thickness:.2f} mm)"
                    )
                    target_modeler.create_box(
                        origin=[
                            "-arrayBoardWidth / 2",
                            "-arrayBoardLength / 2",
                            f"{z_min - z_ref:.4f}mm",
                        ],
                        sizes=["arrayBoardWidth", "arrayBoardLength", f"{thickness:.4f}mm"],
                        name=layer_name,
                        material=material,
                    )

        # --- STEP 4: Replicate Elements via Dedicated Coordinate Systems ---
        replicated_dummy_mapping = {}
        for template in used_templates:
            for key in template.dummies:
                replicated_dummy_mapping.setdefault(key, [])
        replicated_ports = []  # (port name, replicated sheet, CS name)
        element_phase_centres = []
        element_phase_centre_exprs = []
        port_to_cs_mapping = {}

        # Copy each template from its unit cell exactly ONCE. Templates are copied one at a time
        # and tagged '<name>__<template>' right away, so that equal object names in two unit
        # cells (e.g. both have 'L1_Patch') never clash.
        print("\nCopying template geometries to target design...")
        target_modeler.set_working_coordinate_system("Global")
        for template in used_templates:
            names = template.objects
            clashes = [n for n in names if n in target_modeler.object_names]
            if clashes:
                raise RuntimeError(
                    f"Template objects of '{template.design}' clash with objects already in the "
                    f"array design: {', '.join(clashes)}"
                )
            self.target_design_app.copy_solid_bodies_from(
                design=template.app,
                assignment=names,
                no_vacuum=False,
                no_pec=False,
                include_sheets=True,
            )
            missing = [n for n in names if n not in target_modeler.object_names]
            if missing:
                raise RuntimeError(
                    f"Copying templates from '{template.design}' failed for: {', '.join(missing)}"
                )
            for name in names:
                target_modeler.get_object_from_name(name).name = template.tagged(name)
            print(f"  Copied {len(names)} template objects from '{template.design}'.")

        # Clear any copied boundaries/excitations to strip wave ports from template geometries
        try:
            self.target_design_app.oboundary.DeleteAllBoundaries()
            self.target_design_app.oboundary.DeleteAllExcitations()
            print("  Stripped lingering wave port assignments from template sheets.")
        except Exception as e:
            print(f"  Warning: Failed to clear boundaries/excitations ({e})")

        self.target_design_app.set_active_design(self.target_design_name)

        print("Replicating structures to array grid with element coordinate systems...")
        for plan in plans:
            template = plan.template
            label = plan.label
            pos, pos_expr, yaw = plan.position, plan.position_expression, plan.yaw
            print(f"\n  Setting up element {label} ('{template.design}', yaw={yaw:.1f} deg)...")

            # Calculate pointing vectors for the CS rotation around Z
            cos_val, sin_val = math.cos(math.radians(yaw)), math.sin(math.radians(yaw))
            x_pointing = [cos_val, sin_val, 0.0]
            y_pointing = [-sin_val, cos_val, 0.0]

            # Every port gets its own CS at its phase centre (parametric where the layout is),
            # the far-field phase reference of that port. Its origin Z does not affect placement.
            for port in plan.ports:
                if not port.active:
                    continue
                pc = port.phase_centre
                if pos_expr:
                    cs_origin = [
                        offset_expression(pos_expr[0], pc[0] - pos[0]),
                        offset_expression(pos_expr[1], pc[1] - pos[1]),
                        f"{pc[2]:.4f}mm",
                    ]
                else:
                    cs_origin = [f"{pc[0]:.4f}mm", f"{pc[1]:.4f}mm", f"{pc[2]:.4f}mm"]
                element_phase_centres.append(list(pc))
                element_phase_centre_exprs.append(cs_origin[:2])

                # Recreate coordinate system if it exists
                _delete_coordinate_system(target_modeler, port.cs_name)
                print(f"  Creating relative coordinate system '{port.cs_name}' at {cs_origin}...")
                target_modeler.create_coordinate_system(
                    origin=cs_origin,
                    reference_cs="Global",
                    name=port.cs_name,
                    mode="axis",
                    x_pointing=x_pointing,
                    y_pointing=y_pointing,
                )
                # Creating a CS makes it the working CS. All placement operations below are
                # defined in Global coordinates (rotation about the global Z axis, then a global
                # translation), so the working CS must be reset; otherwise the rotation pivots
                # about the element CS and the move vector is read in its rotated axes.
                target_modeler.set_working_coordinate_system("Global")

            # Duplicate the element's template locally using duplicate_along_line with a dummy Z offset of 1.0 mm.
            tagged_to_base = {template.tagged(name): name for name in template.objects}
            success, pasted_names = target_modeler.duplicate_along_line(
                assignment=list(tagged_to_base),
                vector=[0, 0, 1.0],
                clones=2,
            )
            if not success:
                raise RuntimeError(f"Failed to duplicate template objects for element {label}")

            # Move back to origin (offset the Z translation)
            target_modeler.move(
                assignment=pasted_names,
                vector=[0, 0, -1.0],
            )

            # Rename duplicated objects to '<template object>_<label>' and track ports/dummies
            base_names = map_duplicates(pasted_names, tagged_to_base)
            renamed_objs = []
            for pasted_name, base_name in base_names.items():
                new_name = f"{base_name}_{label}"
                target_modeler.get_object_from_name(pasted_name).name = new_name
                renamed_objs.append(new_name)
                for key, dummy_src_list in template.dummies.items():
                    if base_name in dummy_src_list:
                        replicated_dummy_mapping[key].append(new_name)

            # Rotate element if needed (about the global Z axis through the unit-cell origin)
            if yaw != 0.0:
                target_modeler.rotate(
                    assignment=renamed_objs,
                    axis=Axis.Z,
                    angle=yaw,
                )

            # Translate so that the rotated phase centre R*offset lands on the layout position:
            # p -> R(p) + pos - R(offset) = R(p - offset) + pos
            offset = template.offset
            rot_dx = offset[0] * cos_val - offset[1] * sin_val
            rot_dy = offset[0] * sin_val + offset[1] * cos_val

            # In Z, every element is translated by (pos[2] - z_ref), matching the -z_ref shift of
            # the board layers so the stack-up is preserved exactly.
            z_shift_str = f"{pos[2] - z_ref:.4f}mm"

            # Move element to final position
            if pos_expr:
                # Parenthesized so that negative offsets and unary minus in the layout
                # expression never produce ambiguous terms such as "-r - -0.005mm"
                move_vector = [
                    offset_expression(pos_expr[0], -rot_dx),
                    offset_expression(pos_expr[1], -rot_dy),
                    z_shift_str,
                ]
            else:
                move_vector = [
                    f"{pos[0] - rot_dx:.4f}mm",
                    f"{pos[1] - rot_dy:.4f}mm",
                    z_shift_str,
                ]

            target_modeler.move(
                assignment=renamed_objs,
                vector=move_vector,
            )

            # Ports: active ports use their own CS; a matched load uses the element's first
            # active port CS; a dropped port's sheet is removed
            element_cs = next((p.cs_name for p in plan.ports if p.active), "PhaseCentreCS")
            for port in plan.ports:
                if port.name is None:
                    target_modeler.delete([port.sheet])
                    print(f"  Dropped unused port sheet '{port.sheet}'.")
                    continue
                cs = port.cs_name or element_cs
                replicated_ports.append((port.name, port.sheet, cs))
                port_to_cs_mapping[port.name] = cs

        # Reset working coordinate system to Global
        target_modeler.set_working_coordinate_system("Global")
        self.element_phase_centres = element_phase_centres
        self.port_to_cs_mapping = port_to_cs_mapping

        # Clean up the original template geometries at the origin of the target design
        print("\nCleaning up template geometries...")
        target_modeler.delete(
            [template.tagged(name) for template in used_templates for name in template.objects]
        )

        # --- STEP 5: Perform Boolean Dummy Operations ---
        print("\nExecuting boolean cutout operations...")
        for (operation, target_solid), dummy_instances in replicated_dummy_mapping.items():
            if not dummy_instances:
                continue

            # Verify if the target solid exists in the target modeler before executing
            if target_solid not in target_modeler.object_names:
                print(
                    f"  Warning: Target solid '{target_solid}' not found in layout. Skipping boolean {operation}."
                )
                continue

            if operation == "Subtract":
                print(f"  Subtracting {len(dummy_instances)} objects from '{target_solid}'")
                target_modeler.subtract(
                    blank_list=[target_solid],
                    tool_list=dummy_instances,
                    keep_originals=False,
                )
            elif operation == "Unite":
                print(f"  Uniting {len(dummy_instances)} objects with '{target_solid}'")
                target_modeler.unite(assignment=[target_solid] + dummy_instances)

        # --- STEP 6: Assign Wave Port Excitations ---
        print("\nAssigning wave port excitations...")
        # Use first identified ground plane as the default reference ground
        default_ground = ground_layer_names[0] if ground_layer_names else "Ground_Plane"

        # Note on the radiation boundary: the Airbox (STEP 7) is laterally flush with the board, so
        # the wave ports sit on the exterior boundary as HFSS requires (interior wave ports would
        # need a PEC cap). The ABC then terminates the substrate and ground at the board edges,
        # which roughly emulates a laterally extended board for the solve. It is a crude
        # termination: residual reflections limit the reliable XPD floor (~35-40 dB), and the
        # far-field integral treats the dielectric/ground-crossing faces as free space, so
        # near-grazing (theta ~ 90 deg) and back-hemisphere patterns are approximate.
        for port_name, port_sheet, _ in replicated_ports:
            kind = "matched load" if is_load_port(port_name) else "wave port"
            print(f"  Assigning {kind} '{port_name}' to sheet '{port_sheet}' ({default_ground})")
            try:
                self.target_design_app.wave_port(
                    assignment=port_sheet,
                    reference=default_ground,
                    integration_line=Gravity.ZNeg,
                    name=port_name,
                )
            except Exception as e:
                print(
                    f"  Warning: wave_port assignment failed ({e}), falling back to lumped_port..."
                )
                self.target_design_app.lumped_port(
                    assignment=port_sheet,
                    reference=default_ground,
                    integration_line=Gravity.ZNeg,
                    impedance=50.0,
                    name=port_name,
                )

        # Far-field reports show the Edit Sources superposition. Exciting every port at once would
        # mix all Tx/Rx and H/V patterns, so only the first Tx port (or first port) is excited.
        # SBR+ linked antennas use each port separately and are unaffected by this choice.
        port_names = sorted(name for name, _, _ in replicated_ports if not is_load_port(name))
        if port_names:
            tx_ports = [p for p in port_names if "Tx" in p]
            excited_port = tx_ports[0] if tx_ports else port_names[0]
            try:
                # Ports not listed are set to 0 W by PyAEDT
                self.target_design_app.edit_sources(
                    assignment={f"{excited_port}:1": ("1W", "0deg")},
                    include_port_post_processing=True,
                )
                print(
                    f"  Edit Sources: '{excited_port}' at 1 W, all other ports 0 W "
                    "(port post-processing effects included)."
                )
            except Exception as e:
                print(f"  Warning: Could not configure Edit Sources ({e})")

        # --- STEP 6b: Create Array PhaseCentreCS & Configure Far-field Phase Reference ---
        # PhaseCentreCS sits at the centroid of the port phase centres with global axes. It is
        # the CS of all far-field spheres (z = boresight, x = co-polar reference). The far-field
        # phase reference for SBR+ linking is per port (port CSs) unless mode 'single' puts it
        # here. Lateral coordinates stay parametric when the layout is.
        if element_phase_centres:
            count = len(element_phase_centres)
            avg_z = sum(pt[2] for pt in element_phase_centres) / count
            pc_origin = [
                f"({' + '.join(e[axis] for e in element_phase_centre_exprs)}) / {count}"
                for axis in (0, 1)
            ] + [f"{avg_z:.4f}mm"]
        else:
            pc_origin = ["0mm", "0mm", "0mm"]

        _delete_coordinate_system(target_modeler, "PhaseCentreCS")

        print(f"\nCreating array PhaseCentreCS at element centroid {pc_origin}...")
        target_modeler.create_coordinate_system(
            origin=pc_origin,
            reference_cs="Global",
            name="PhaseCentreCS",
            mode="axis",
            x_pointing=[1.0, 0.0, 0.0],
            y_pointing=[0.0, 1.0, 0.0],
        )
        target_modeler.set_working_coordinate_system("Global")

        # Set Far-field Phase Reference for SBR+ Antenna Linking
        self.set_far_field_phase_reference(mode=self.phase_reference_mode)

        # --- STEP 7: Create Radiation Airbox ---
        print("\nCreating radiation boundary Airbox...")
        # Get overall board bounds in Z direction
        all_z_coords = []
        for layer_info in global_layers.values():
            all_z_coords.extend([layer_info["z_min"], layer_info["z_max"]])
        overall_z_min = (min(all_z_coords) if all_z_coords else 0.0) - z_ref
        overall_z_max = (max(all_z_coords) if all_z_coords else 1.0) - z_ref
        total_z_thickness = overall_z_max - overall_z_min

        print(
            f"  Creating flush lateral Airbox solid (airboxClearance={airbox_clearance_mm:.2f} mm)"
        )
        airbox_obj = target_modeler.create_box(
            origin=[
                "-arrayBoardWidth / 2",
                "-arrayBoardLength / 2",
                f"{overall_z_min:.4f}mm - airboxClearance",
            ],
            sizes=[
                "arrayBoardWidth",
                "arrayBoardLength",
                f"{total_z_thickness:.4f}mm + 2 * airboxClearance",
            ],
            name="Airbox",
            material="vacuum",
        )
        airbox_obj.display_wireframe = True

        print("  Assigning radiation boundary to Airbox...")
        try:
            self.target_design_app.assign_radiation_boundary_to_objects("Airbox", "Radiation_Box")
        except Exception as e:
            print(f"  Warning: Failed to assign radiation boundary ({e})")
        try:
            airbox_obj.transparency = 1.0
        except Exception as e:
            print(f"  Warning: Failed to set Airbox transparency ({e})")

        # --- STEP 7b: Post-build placement check ---
        placement_errors = self.verify_placement(
            plans, z_ref, [sheet for _, sheet, _ in replicated_ports]
        )
        if placement_errors:
            self.target_design_app.save_project()
            raise RuntimeError(
                "Post-build placement check failed (project saved for inspection):\n  "
                + "\n  ".join(placement_errors)
            )

        # --- STEP 7c: Port polarization report ---
        self.report_port_polarizations(plans, z_ref)

        # Configure simulation setup and sweep by inheriting from the element design
        self.configure_simulation_setup()

        # Far-field spheres and post-processing reports
        self.create_post_processing_reports()

        # Run Design Validation
        print("\nRunning HFSS built-in design validation...")
        validation_ok = self.target_design_app.validate_simple()
        if validation_ok == 1 or validation_ok is True:
            print("  Design validation: PASSED.")
        else:
            self.target_design_app.save_project()
            raise RuntimeError(
                "HFSS design validation failed (project saved for inspection in the GUI)."
            )

        self.target_design_app.save_project()
        print(
            f"\nMIMO array synthesis successfully completed in design '{self.target_design_name}'."
        )

        self.run_solve(simulate=simulate)

    @staticmethod
    def _template_corners(template):
        """
        Footprint corners of a template in its own coordinates: its layer outline in x (widened
        to any copper reaching beyond it, e.g. the side feed of a dual-port cell) and the reach of
        its active structures (feedlines and ports) in y.
        """
        modeler = template.app.modeler
        x_values, y_layer, y_active = [], [], []
        for name in template.layers:
            bbox = modeler.get_object_from_name(name).bounding_box
            x_values += [float(bbox[0]), float(bbox[3])]
            y_layer += [float(bbox[1]), float(bbox[4])]
        for name in template.active + template.port_sheets:
            try:
                bbox = modeler.get_object_from_name(name).bounding_box
                y_active += [float(bbox[1]), float(bbox[4])]
                if x_values:
                    x_values += [float(bbox[0]), float(bbox[3])]
            except Exception:
                pass
        y_values = y_active or y_layer
        if not x_values or not y_values:
            return []
        return [(x, y) for x in (min(x_values), max(x_values)) for y in (min(y_values), max(y_values))]

    def report_port_polarizations(self, plans, z_ref):
        """
        Prints every port's built polarization and field sign, warns about sign combinations
        that a sphere calibration cannot reveal, and writes '<design>_ports.json' next to the
        project.
        """
        rows = port_table(plans)
        print("\nPort polarizations (global frame: H along +x, V along +y):")
        print(f"  {'Port':<24} {'Role':<4} {'Pol':<3} {'Sign':>4}  {'Template':<24} {'Yaw':>6}")
        for row in rows:
            port = row["port"] or f"({row['template_sheet']} dropped)"
            sign = {1: "+", -1: "-"}.get(row["sign"], "?")
            kind = "" if row["active"] else "  load"
            print(
                f"  {port:<24} {row['role'] or '?':<4} {row['polarization'] or '?':<3} "
                f"{sign:>4}  {row['template']:<24} {row['yaw_deg']:>6g}{kind}"
            )
        for warning in sign_warnings(rows):
            print(f"  [WARNING] {warning}")

        sidecar = os.path.join(
            os.path.dirname(self.project_path), f"{self.target_design_name}_ports.json"
        )
        try:
            with open(sidecar, "w", encoding="utf-8") as f:
                json.dump(sidecar_dict(self.target_design_name, plans, z_ref), f, indent=2)
            print(f"  Port table written to {sidecar}")
        except OSError as e:
            print(f"  Warning: Could not write the port table ({e})")

    def verify_placement(self, plans, z_ref, port_sheet_names, tolerance_mm=1e-3):
        """
        Checks the built geometry against the layout, independently of how AEDT applied the
        placement operations.

        1. Every replicated template object must satisfy p -> R(yaw) (p - offset_xy) + pos_xy and
           z -> z - z_ref + pos_z (its template's offset, the common stack reference), compared
           via the vertex centroid (exact under any rotation).
        2. Every wave port sheet must lie on a lateral face of the Airbox (exterior boundary).

        Returns a list of human-readable violations (empty when everything is in place).
        """
        target_modeler = self.target_design_app.modeler
        print("\nVerifying element placement and port positions...")

        def vertex_centroid(obj):
            positions = [vertex.position for vertex in obj.vertices]
            return np.mean(np.array(positions, dtype=float), axis=0) if positions else None

        errors = []
        centroids = {}
        for plan in plans:
            template = plan.template
            if id(template) in centroids:
                continue
            template_centroids = {}
            for name in template.active + template.port_sheets:
                centroid = vertex_centroid(template.app.modeler[name])
                if centroid is None:
                    print(f"  Note: '{name}' has no vertices; its placement is not checked.")
                else:
                    template_centroids[name] = centroid
            centroids[id(template)] = template_centroids

        checked = 0
        for plan in plans:
            dropped = {port.spec.sheet for port in plan.ports if port.name is None}
            for name, template_centroid in centroids[id(plan.template)].items():
                if name in dropped:
                    continue
                obj = target_modeler.get_object_from_name(f"{name}_{plan.label}")
                if not obj:
                    errors.append(f"Missing replicated object '{name}_{plan.label}'")
                    continue
                expected = np.array(
                    place(template_centroid, plan.template.offset, z_ref, plan.yaw, plan.position)
                )
                actual = vertex_centroid(obj)
                deviation = np.linalg.norm(actual - expected)
                checked += 1
                if deviation > tolerance_mm:
                    errors.append(
                        f"'{name}_{plan.label}' is {deviation:.4f} mm off its layout position "
                        f"(expected centroid {np.round(expected, 4).tolist()}, "
                        f"actual {np.round(actual, 4).tolist()})"
                    )

        airbox = target_modeler.get_object_from_name("Airbox")
        if airbox:
            airbox_bbox = [float(v) for v in airbox.bounding_box]
            for sheet_name in port_sheet_names:
                bbox = [float(v) for v in target_modeler[sheet_name].bounding_box]
                flat_axes = [
                    axis for axis in (0, 1) if abs(bbox[axis] - bbox[axis + 3]) < tolerance_mm
                ]
                on_face = any(
                    min(
                        abs(bbox[axis] - airbox_bbox[axis]),
                        abs(bbox[axis] - airbox_bbox[axis + 3]),
                    )
                    < tolerance_mm
                    for axis in flat_axes
                )
                if not on_face:
                    errors.append(
                        f"Port sheet '{sheet_name}' does not lie on a lateral Airbox face "
                        f"(bbox {np.round(bbox, 4).tolist()}); interior wave ports need a PEC cap"
                    )

        if errors:
            print(f"  [ERROR] {len(errors)} placement violation(s) found.")
        else:
            print(
                f"  OK: {checked} replicated objects within {tolerance_mm * 1000:.0f} um of their "
                f"layout positions; {len(port_sheet_names)} port sheets on the Airbox boundary."
            )
        return errors

    def configure_simulation_setup(self):
        """
        Gives the array design 'Setup': the primary unit cell's analysis setup for its
        meshing (adaptive frequencies, convergence and solver settings, as validated on the unit
        cell), with the copied sweeps replaced by the builder's own discrete sweep over the
        frequencies of the unit cell's sweep (`_create_array_sweep`), which SBR+ can link. The
        setup is copy-pasted natively, else recreated from the unit cell's properties, else the
        build fails; a unit cell without a sweep fails before anything is changed.
        """
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

        # Clear any leftover Optimetrics setups that may be unlinked
        try:
            opt_names = list(self.target_design_app.ooptimetrics.GetSetupNames())
            if opt_names:
                print(f"  Clearing {len(opt_names)} existing Optimetrics setups...")
                self.target_design_app.ooptimetrics.DeleteSetups(opt_names)
        except Exception:
            pass

        # Check and purge any existing analysis setups and their sweeps
        try:
            setup_names = list(self.target_design_app.setup_names)
        except Exception:
            setup_names = []

        if setup_names:
            print(f"  Clearing {len(setup_names)} existing analysis setups and their sweeps in target design...")
            for s_name in setup_names:
                try:
                    setup_obj = self.target_design_app.get_setup(s_name)
                    # Initialize PyAEDT sweeps cache to prevent NoneType iterable error
                    try:
                        _ = setup_obj.sweeps
                    except Exception:
                        pass
                    # Delete all sweeps inside this setup
                    try:
                        sweep_names = list(setup_obj.get_sweep_names())
                    except Exception:
                        sweep_names = []
                    for sw_name in sweep_names:
                        try:
                            setup_obj.delete_sweep(sw_name)
                            print(f"  Deleted existing sweep '{sw_name}' from setup '{s_name}'")
                        except Exception:
                            pass
                except Exception:
                    pass
                try:
                    self.target_design_app.delete_setup(s_name)
                    print(f"  Deleted existing setup '{s_name}'")
                except Exception:
                    pass

        target_setup_name = ARRAY_SETUP_NAME
        source_setups = list(self.source_design_app.setup_names) if self.source_design_app else []
        if not source_setups:
            raise RuntimeError(
                f"Unit cell '{self.source_design_name}' has no analysis setup to copy; the array "
                "takes its meshing and frequency sweep from it."
            )
        src_setup_name = getattr(self.source_design_app, "active_setup", None) or source_setups[0]
        src_setup = self.source_design_app.get_setup(src_setup_name)
        # The array's frequencies come from the unit cell's (first) sweep, as saved
        saved = saved_setups(self.source_design_app).get(src_setup_name, {})
        src_sweeps = list(self.source_design_app.oanalysis.GetSweeps(src_setup_name))
        if not src_sweeps:
            raise RuntimeError(
                f"Unit cell setup '{src_setup_name}' has no frequency sweep; the array's sweep "
                "takes its frequencies from it."
            )
        if len(src_sweeps) > 1:
            print(f"  Note: the unit cell has sweeps {src_sweeps}; the first sets the frequencies.")
        origin = f"{src_setup_name} : {src_sweeps[0]}"
        ranges = sweep_ranges(setup_sweeps(saved).get(src_sweeps[0], {}))
        if not ranges:
            raise RuntimeError(f"Cannot read the frequency range of the unit cell's '{origin}'.")

        print(f"\nCopying the unit cell's analysis setup '{src_setup_name}' for its meshing...")
        self._copy_unit_cell_setup(src_setup_name, src_setup, target_setup_name)
        self._create_array_sweep(target_setup_name, ranges, origin)
        coverage = adaptive_coverage_warning(saved, ranges)
        if coverage:
            print(f"  Warning: {coverage}")

    def _copy_unit_cell_setup(self, src_setup_name, src_setup, target_setup_name):
        """
        Copies the unit cell's setup into the array design as `target_setup_name`: natively
        (copy-paste, which brings its sweeps along; `_create_array_sweep` replaces them), else
        recreated from its adaptive and solver properties, else raises.
        """
        try:
            self.source_design_app.oanalysis.CopyDrivenSetup(src_setup_name)
            pasted_name = self.target_design_app.oanalysis.PasteDrivenSetup()
            if pasted_name:
                if pasted_name != target_setup_name:
                    self.target_design_app.oanalysis.RenameSetup(pasted_name, target_setup_name)
                print(f"  Copy-pasted the setup as '{target_setup_name}'.")
                return
        except Exception as err:
            print(
                f"  Warning: Native CopyDrivenSetup/PasteDrivenSetup encountered an issue ({err}). Falling back to manual recreation..."
            )

        try:
            print(f"  Recreating setup '{target_setup_name}' from source properties...")
            solve_type = str(src_setup.props.get("SolveType", "BroadBand")).lower()

            setup = self.target_design_app.create_setup(name=target_setup_name)
            setup.auto_update = False

            max_passes = src_setup.props.get("MaximumPasses", 21)
            max_delta_s = src_setup.props.get("MaxDeltaS", 0.02)

            if solve_type == "broadband":
                adaptive_freqs = src_setup.props.get("MultipleAdaptiveFreqsSetup", {})
                low_f, high_f = adaptive_freqs["Low"], adaptive_freqs["High"]
                setup.enable_adaptive_setup_broadband(
                    low_frequency=low_f,
                    high_frequency=high_f,
                    max_passes=max_passes,
                    max_delta_s=max_delta_s,
                )
                print(f"  Configured Broadband adaptive mesh: Low={low_f}, High={high_f}, MaxPasses={max_passes}, MaxDeltaS={max_delta_s}")
            elif solve_type == "multifrequency":
                freq_list = [f"{f:g}GHz" for f in adaptive_frequencies_ghz(src_setup.props)]
                if not freq_list:
                    raise RuntimeError("cannot read the multi-frequency adaptive frequencies")
                setup.enable_adaptive_setup_multifrequency(
                    frequencies=freq_list,
                    max_delta_s=max_delta_s,
                )
                print(f"  Configured Multi-Frequency adaptive mesh: Frequencies={freq_list}, MaxDeltaS={max_delta_s}")
            else:
                freq = src_setup.props.get("Frequency", f"{self.centre_frequency_ghz}GHz")
                setup.enable_adaptive_setup_single(
                    freq=freq,
                    max_passes=max_passes,
                    max_delta_s=max_delta_s,
                )
                print(f"  Configured Single-Frequency adaptive mesh: Frequency={freq}, MaxPasses={max_passes}, MaxDeltaS={max_delta_s}")

            # Replicate meshing, convergence, and solver options
            for prop_key in [
                "BasisOrder", "DoLambdaRefine", "DoMaterialLambda", "SetLambdaTarget",
                "Target", "UseMaxTetIncrease", "PortAccuracy", "SaveFields", "SaveAnyFields",
                "PercentRefinement", "MinimumPasses", "MinimumConvergedPasses", "DrivenSolverType"
            ]:
                if prop_key in src_setup.props:
                    setup.props[prop_key] = src_setup.props[prop_key]

            setup.auto_update = True
            setup.update()
        except Exception as err2:
            raise RuntimeError(
                f"Could not copy the analysis setup '{src_setup_name}' of unit cell "
                f"'{self.source_design_name}' ({err2})."
            ) from err2

    def _create_array_sweep(self, setup_name, ranges, origin):
        """
        Replaces every sweep of the array setup `setup_name` by the builder's own: one discrete
        sweep 'Sweep' over the frequency `ranges` of the unit cell's sweep `origin`, saving the
        radiated fields at every frequency. SBR+ receives a linked antenna's far fields only at
        frequencies where it saved fields and interpolates in between, so an interpolating
        sweep (fields at its basis frequencies) or the adaptive solution would put ripple into
        every SBR+ result (see frequency_sweeps). A discrete sweep solves the whole array once
        per frequency on the final adaptive mesh.
        """
        analysis = self.target_design_app.oanalysis
        for name in list(analysis.GetSweeps(setup_name)):
            analysis.DeleteSweep(setup_name, name)
            print(f"  Removed the copied sweep '{name}'.")
        settings = [
            ("IsEnabled", True),
            ("Type", LINKABLE_SWEEP_TYPE),
            ("SaveFields", True),
            ("SaveRadFields", True),
            ("ExtrapToDC", False),
        ]
        analysis.InsertFrequencySweep(
            setup_name, native_sweep_arguments(ARRAY_SWEEP_NAME, ranges, settings)
        )
        if ARRAY_SWEEP_NAME not in list(analysis.GetSweeps(setup_name)):
            raise RuntimeError(f"AEDT did not create '{setup_name} : {ARRAY_SWEEP_NAME}'.")
        low, high = sweep_bounds_ghz(ranges)
        try:
            points = f"{sweep_point_count(ranges)} frequencies"
        except ValueError:
            points = "its frequencies"
        print(
            f"  Created the discrete sweep '{setup_name} : {ARRAY_SWEEP_NAME}' over "
            f"{low:g}-{high:g} GHz ({points}, from the unit cell's '{origin}'), saving the "
            "radiated fields; solving it runs one full-wave solve of the array per frequency."
        )

    def create_post_processing_reports(self):
        """
        Creates standard S-parameter and Far-Field reports in the target design.

        Existing S-parameter reports are preserved. The builder-owned far-field spheres and
        reports are enforced/recreated (see `configure_far_field_spheres`): realized gain is used
        for absolute patterns and directivity for XPD.
        """
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

        if self.centre_frequency_ghz is None:
            self.centre_frequency_ghz = (
                _setup_frequency_ghz(self.target_design_app) or DEFAULT_CENTRE_FREQUENCY_GHZ
            )

        print("\nCreating automated post-processing reports...")
        try:
            raw_excitations = self.target_design_app.excitation_names
            # Strip mode suffixes (e.g. "Port_Rx_1_V:1" -> "Port_Rx_1_V") and deduplicate;
            # matched loads on unused feeds are left out of the reports
            port_names = sorted(
                set(p.split(":")[0] for p in raw_excitations if p and not is_load_port(p))
            )
        except Exception as e:
            print(f"  Warning: Failed to retrieve excitation names ({e})")
            port_names = []

        if not port_names:
            print("  Warning: No ports found. Skipping S-parameter reports.")
            return

        rx_ports = sorted([p for p in port_names if "Rx" in p])
        tx_ports = sorted([p for p in port_names if "Tx" in p])

        target_setups = list(self.target_design_app.setup_names)
        if target_setups:
            s_name = ARRAY_SETUP_NAME if ARRAY_SETUP_NAME in target_setups else target_setups[0]
            try:
                s_obj = self.target_design_app.get_setup(s_name)
                sw_names = s_obj.get_sweep_names()
                sw_name = sw_names[0] if sw_names else "Sweep"
            except Exception:
                sw_name = "Sweep"
            setup_sweep = f"{s_name} : {sw_name}"
        else:
            setup_sweep = f"{ARRAY_SETUP_NAME} : {ARRAY_SWEEP_NAME}"

        # Get list of existing reports to avoid duplicate report generation
        try:
            existing_reports = list(self.target_design_app.post.all_report_names)
        except Exception:
            existing_reports = []

        # --- STEP A: Set up S-parameters FIRST ---
        # 1. Reflections S-parameters (S_ii)
        plot_name = "Reflections"
        if plot_name not in existing_reports:
            print("  Generating Reflections S-parameter report...")
            reflections = [f"dB(S({p},{p}))" for p in port_names]
            try:
                self.target_design_app.post.create_report(
                    expressions=reflections,
                    setup_sweep_name=setup_sweep,
                    plot_name=plot_name,
                    report_category="Modal Solution Data",
                )
            except Exception as e:
                print(f"  Warning: Failed to create Reflections report ({e})")
        else:
            print(f"  Preserving existing report '{plot_name}'")

        # 2. Rx Crosstalk
        if len(rx_ports) > 1:
            plot_name = "Rx crosstalk"
            if plot_name not in existing_reports:
                print("  Generating Rx crosstalk report...")
                rx_couplings = []
                for i in range(len(rx_ports)):
                    for j in range(i):
                        rx_couplings.append(f"dB(S({rx_ports[i]},{rx_ports[j]}))")
                try:
                    self.target_design_app.post.create_report(
                        expressions=rx_couplings,
                        setup_sweep_name=setup_sweep,
                        plot_name=plot_name,
                        report_category="Modal Solution Data",
                    )
                except Exception as e:
                    print(f"  Warning: Failed to create Rx crosstalk report ({e})")
            else:
                print(f"  Preserving existing report '{plot_name}'")

        # 3. Tx Crosstalk
        if len(tx_ports) > 1:
            plot_name = "Tx crosstalk"
            if plot_name not in existing_reports:
                print("  Generating Tx crosstalk report...")
                tx_couplings = []
                for i in range(len(tx_ports)):
                    for j in range(i):
                        tx_couplings.append(f"dB(S({tx_ports[i]},{tx_ports[j]}))")
                try:
                    self.target_design_app.post.create_report(
                        expressions=tx_couplings,
                        setup_sweep_name=setup_sweep,
                        plot_name=plot_name,
                        report_category="Modal Solution Data",
                    )
                except Exception as e:
                    print(f"  Warning: Failed to create Tx crosstalk report ({e})")
            else:
                print(f"  Preserving existing report '{plot_name}'")

        # 4. Tx1-to-Rx Crosstalk
        if tx_ports and rx_ports:
            plot_name = "Tx1-to-Rx crosstalk"
            if plot_name not in existing_reports:
                tx1 = tx_ports[0]
                print("  Generating Tx1-to-Rx crosstalk report...")
                tx_to_rx = [f"dB(S({rx},{tx1}))" for rx in rx_ports]
                try:
                    self.target_design_app.post.create_report(
                        expressions=tx_to_rx,
                        setup_sweep_name=setup_sweep,
                        plot_name=plot_name,
                        report_category="Modal Solution Data",
                    )
                except Exception as e:
                    print(f"  Warning: Failed to create Tx1-to-Rx crosstalk report ({e})")
            else:
                print(f"  Preserving existing report '{plot_name}'")

        # --- STEP B: Far-field spheres and radiation pattern reports ---
        # Ensure far-field phase reference is configured for SBR+ linking
        try:
            self.set_far_field_phase_reference(mode=self.phase_reference_mode)
        except Exception as e:
            print(f"  Warning: Failed to configure far-field phase reference ({e})")

        # Builder-owned far-field reports are recreated, since their sweep variables depend on
        # the sphere definitions (enforced below). Reports with other names are left untouched.
        for plot_name in (
            FF_REPORT_3D,
            FF_REPORT_PHI0,
            FF_REPORT_PHI90,
            FF_REPORT_XPD_PHI0,
            FF_REPORT_XPD_PHI90,
        ):
            if plot_name in existing_reports:
                self.target_design_app.post.delete_report(plot_name)

        self.configure_far_field_spheres()

        freq = f"{self.centre_frequency_ghz}GHz"
        cut_variations = {"Theta": ["All"], "Phi": ["All"], "Freq": [freq]}
        report_specs = [
            # (name, sphere, expressions, extra create_report kwargs)
            (
                FF_REPORT_3D,
                FF_SPHERE_3D,
                ["dB(RealizedGainTotal)"],
                dict(
                    primary_sweep_variable="Phi",
                    secondary_sweep_variable="Theta",
                    plot_type="3D Polar Plot",
                ),
            ),
            (
                FF_REPORT_PHI0,
                FF_SPHERE_PHI0,
                ["dB(RealizedGainTotal)", "dB(RealizedGainCoPolar)", "dB(RealizedGainCrossPolar)"],
                dict(primary_sweep_variable="Theta", plot_type="Rectangular Plot"),
            ),
            (
                FF_REPORT_PHI90,
                FF_SPHERE_PHI90,
                ["dB(RealizedGainTotal)", "dB(RealizedGainCoPolar)", "dB(RealizedGainCrossPolar)"],
                dict(primary_sweep_variable="Theta", plot_type="Rectangular Plot"),
            ),
            (
                FF_REPORT_XPD_PHI0,
                FF_SPHERE_PHI0,
                [XPD_X_POL_EXPRESSION, XPD_Y_POL_EXPRESSION],
                dict(primary_sweep_variable="Theta", plot_type="Rectangular Plot"),
            ),
            (
                FF_REPORT_XPD_PHI90,
                FF_SPHERE_PHI90,
                [XPD_X_POL_EXPRESSION, XPD_Y_POL_EXPRESSION],
                dict(primary_sweep_variable="Theta", plot_type="Rectangular Plot"),
            ),
        ]
        for plot_name, sphere_name, expressions, kwargs in report_specs:
            print(f"  Generating '{plot_name}' report on '{sphere_name}' at {freq}...")
            try:
                self.target_design_app.post.create_report(
                    expressions=expressions,
                    setup_sweep_name=setup_sweep,
                    variations=dict(cut_variations),
                    report_category="Far Fields",
                    plot_name=plot_name,
                    context=sphere_name,
                    **kwargs,
                )
            except Exception as e:
                print(f"  Warning: Failed to create '{plot_name}' report ({e})")

    def configure_far_field_spheres(self, cs_name="PhaseCentreCS"):
        """
        Creates the builder-owned far-field spheres in the array design, or overwrites their
        definition in place when they exist (so that user reports on them survive).

        All spheres use the Theta-Phi definition in `cs_name` with an explicit Z-axis boresight and
        Linear polarization. CoPolar/CrossPolar are then Ludwig-3 X/Y referenced to the CS x-axis
        on every sphere. Ludwig-3 is regular at boresight (theta = 0) and singular only at
        theta = 180 deg. The principal cuts use signed theta in [-180, 180] deg at fixed phi:
        positive theta lies in the +x (phi = 0) or +y (phi = 90) half-plane.
        """
        app = self.target_design_app
        cs_names = {cs.name for cs in app.modeler.coordinate_systems}
        sphere_cs = cs_name if cs_name in cs_names else None
        if sphere_cs is None:
            print(f"  Warning: '{cs_name}' not found; far-field spheres use the Global CS.")

        grids = {
            FF_SPHERE_3D: ((0, 180, 1), (0, 360, 5)),
            FF_SPHERE_PHI0: ((-180, 180, 1), (0, 0, 1)),
            FF_SPHERE_PHI90: ((-180, 180, 1), (90, 90, 1)),
        }
        existing = {setup.name for setup in app.field_setups}
        for name, ((t_start, t_stop, t_step), (p_start, p_stop, p_step)) in grids.items():
            props = {
                "UseCustomRadiationSurface": False,
                "CSDefinition": "Theta-Phi",
                "Polarization": "Linear",
                "SlantAngle": "0deg",
                "Boresight": "Z Axis",
                "ThetaStart": f"{t_start}deg",
                "ThetaStop": f"{t_stop}deg",
                "ThetaStep": f"{t_step}deg",
                "PhiStart": f"{p_start}deg",
                "PhiStop": f"{p_stop}deg",
                "PhiStep": f"{p_step}deg",
                "UseLocalCS": sphere_cs is not None,
                "CoordSystem": sphere_cs or "",
            }
            args = ["NAME:" + name]
            _dict2arg(props, args)
            try:
                if name in existing:
                    app.oradfield.EditInfiniteSphereSetup(name, args)
                    action = "Updated"
                else:
                    app.oradfield.InsertInfiniteSphereSetup(args)
                    action = "Created"
                print(
                    f"  {action} far-field sphere '{name}' "
                    f"(Theta-Phi in '{sphere_cs or 'Global'}', "
                    f"theta {t_start}..{t_stop} deg, phi {p_start}..{p_stop} deg, boresight Z)"
                )
            except Exception as e:
                print(f"  Warning: Failed to configure far-field sphere '{name}' ({e})")

    def run_solve(self, simulate=None):
        """Solves the array setup (adaptive passes and the discrete sweep), after asking."""
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

        question = (
            f"Start the full-wave simulation of '{self.target_design_name}' now? The discrete "
            "sweep solves the whole array once per frequency"
        )
        if not _confirm(question, preset=simulate, default=False):
            print("Simulation not started.")
            return
        # The setup is deliberately not reconfigured here, since that deletes existing setups
        # together with their solutions
        print(f"Solving HFSS array design '{self.target_design_name}'...")
        self.target_design_app.analyze_setup(ARRAY_SETUP_NAME)
        self.target_design_app.save_project()

    def export_coupling_s_parameters(self, output_touchstone_path, setup_name=ARRAY_SETUP_NAME):
        """Exports solved S-parameter coupling matrix to a Touchstone (.sNp) file."""
        if not self.target_design_app:
            raise RuntimeError("Array design is not resolved.")

        output_touchstone_path = os.path.abspath(output_touchstone_path)
        sweep_names = self.target_design_app.get_setup(setup_name).get_sweep_names()
        sweep_name = sweep_names[0] if sweep_names else "LastAdaptive"
        print(
            f"Exporting coupling S-parameters ({setup_name} : {sweep_name}) "
            f"to: {output_touchstone_path}"
        )

        self.target_design_app.export_touchstone(
            setup=setup_name, sweep=sweep_name, output_file=output_touchstone_path
        )

    def set_far_field_phase_reference(self, mode=None, cs_name="PhaseCentreCS"):
        """
        Sets the Far-Field Phase Reference in the array design
        (HFSS -> Excitations -> Set Far-field Phase Reference).
        Configures the phase reference origin used when linking the antenna to SBR+ simulations.

        Parameters
        ----------
        mode : str, optional
            Phase reference mode:
            - "per_port" or "per-port" (default): each port's far field is referenced to its own
              element CS 'CS_<label>' (at the element phase centre). A linked SBR+ antenna then
              launches each port's rays from that point, so the Tx/Rx geometry (bistatic angles,
              path lengths) is modelled by the ray tracer and stays exact for near targets.
            - "single": all ports share 'PhaseCentreCS' (array centroid). SBR+ then launches every
              port from the centroid and the element offsets survive only as far-field pattern
              phase, which is valid only for targets in the far field of the whole array.
        cs_name : str, optional
            The coordinate system name to use when mode="single". Default is "PhaseCentreCS".
        """
        if not self.target_design_app:
            raise RuntimeError("Target array design not connected or synthesized yet.")

        if mode is None:
            mode = getattr(self, "phase_reference_mode", "per_port")
        normalized_mode = str(mode).lower().replace("-", "_")

        print(f"\nConfiguring Far-field Phase Reference (mode='{normalized_mode}')...")

        try:
            cs_names = {cs.name for cs in self.target_design_app.modeler.coordinate_systems}
        except Exception:
            cs_names = set()

        if normalized_mode == "single":
            target_cs = cs_name if cs_name in cs_names else "Global"
            try:
                self.target_design_app.oboundary.SetSinglePhaseCenter(target_cs)
                print(f"  Successfully set Single Far-field Phase Reference to '{target_cs}'.")
            except Exception as e:
                print(f"  Warning: Failed to set single phase centre to '{target_cs}' ({e})")
        elif normalized_mode == "per_port":
            try:
                raw_excitations = self.target_design_app.excitation_names
                port_names = sorted(list(set(p.split(":")[0] for p in raw_excitations if p)))
            except Exception:
                port_names = []

            mapping = []
            for port in port_names:
                assigned_cs = None
                if hasattr(self, "port_to_cs_mapping") and port in self.port_to_cs_mapping:
                    assigned_cs = self.port_to_cs_mapping[port]
                elif port.startswith("Port_") and f"CS_{port[5:]}" in cs_names:
                    assigned_cs = f"CS_{port[5:]}"
                elif f"CS_{port}" in cs_names:
                    assigned_cs = f"CS_{port}"
                elif cs_name in cs_names:
                    assigned_cs = cs_name
                else:
                    assigned_cs = "Global"

                mapping.append(["NAME:" + port, "Coordinate System:=", assigned_cs])

            try:
                self.target_design_app.oboundary.SetPhaseCenterPerPort(mapping)
                summary_str = ", ".join(
                    f"{m[0].replace('NAME:', '')} -> {m[2]}" for m in mapping
                )
                print(
                    f"  Successfully set Per-Port Far-field Phase Reference for {len(mapping)} ports: {summary_str}"
                )
            except Exception as e:
                print(f"  Warning: Failed to set per-port phase centre ({e})")
        else:
            print(f"  Warning: Unknown phase reference mode '{mode}'. Skipping.")

    def fit_unit_cell_phase_centre(self, cap_deg=40.0, cs_name="PhaseCentreCS", app=None):
        """
        Fits the unit-cell phase centre from its solved far field and stores it as `cs_name` in the
        unit-cell design.

        The co-polar Ludwig-3 component (L3X or L3Y, whichever dominates at boresight) is sampled
        over the cap theta <= cap_deg (all phi) at the centre frequency. A phase centre at d makes
        the far-field phase vary as arg E_co(u) = psi_0 + k u.d (u = unit direction), so d is the
        |E_co|-weighted least-squares solution of that linear model. Unlike a single-cut phase
        flatness optimization, this constrains X, Y and Z and needs no extra solves. The separate
        phi = 0 / 90 plane fits are reported as a measure of astigmatism.

        Returns
        -------
        list of float
            Phase centre [x, y, z] in mm, in unit-cell global coordinates.
        """
        app = app or self.source_design_app
        freq_ghz = (
            self.centre_frequency_ghz or _setup_frequency_ghz(app) or DEFAULT_CENTRE_FREQUENCY_GHZ
        )
        print(
            f"\nFitting unit-cell phase centre (theta <= {cap_deg:.0f} deg at {freq_ghz:.3f} GHz)..."
        )

        sphere_name = "PhaseCentreFitSphere"
        _delete_far_field_sphere(app, sphere_name)
        app.insert_infinite_sphere(
            name=sphere_name,
            definition="Theta-Phi",
            theta_start=0,
            theta_stop=cap_deg,
            theta_step=2,
            phi_start=0,
            phi_stop=355,
            phi_step=5,
            units="deg",
        )
        try:
            data = self._unit_cell_far_field(sphere_name, ["rEL3X", "rEL3Y"], freq_ghz, app)
        finally:
            _delete_far_field_sphere(app, sphere_name)

        theta_deg, phi_deg = data["Theta"], data["Phi"]
        boresight = np.isclose(theta_deg, 0.0)
        co_name = max(("rEL3X", "rEL3Y"), key=lambda e: np.mean(np.abs(data[e][boresight])))
        e_co = data[co_name]

        theta, phi = np.radians(theta_deg), np.radians(phi_deg)
        u = np.column_stack(
            [np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]
        )
        k = 2.0 * np.pi * freq_ghz * 1e9 / SPEED_OF_LIGHT_MM_S  # rad/mm
        weights = np.abs(e_co)

        # Phase relative to boresight, unwrapped outward along theta on every phi cut
        psi = np.angle(e_co / np.mean(e_co[boresight]))
        for phi_value in np.unique(phi_deg):
            idx = np.flatnonzero(np.isclose(phi_deg, phi_value))
            idx = idx[np.argsort(theta_deg[idx])]
            psi[idx] = np.unwrap(psi[idx])

        def solve(mask, columns):
            design = np.column_stack([np.ones(mask.sum()), k * u[mask][:, columns]])
            solution, *_ = np.linalg.lstsq(
                design * weights[mask, None], psi[mask] * weights[mask], rcond=None
            )
            return solution[1:]

        def rms_phase_error_deg(centre):
            shifted = e_co * np.exp(-1j * k * (u @ centre))
            mean = np.sum(weights * shifted) / np.sum(weights)
            error = np.angle(shifted / mean)
            return np.degrees(np.sqrt(np.sum(weights**2 * error**2) / np.sum(weights**2)))

        everywhere = np.ones_like(psi, dtype=bool)
        centre = solve(everywhere, [0, 1, 2])
        phi0_cut = np.isclose(np.mod(phi_deg, 180.0), 0.0)
        phi90_cut = np.isclose(np.mod(phi_deg, 180.0), 90.0)
        z_phi0 = solve(phi0_cut, [0, 2])[1]
        z_phi90 = solve(phi90_cut, [1, 2])[1]

        centre_mm = [round(float(c), 3) for c in centre]
        print(f"  Co-polar component: {co_name}")
        print(f"  Fitted phase centre: {centre_mm} mm")

        # Lateral offsets shift each (rotated) element against the symmetric board outline, which
        # would pull the wave ports off the flush radiation boundary. Sub-10 um lateral components
        # are fit noise (< 1 deg phase over the cap at mm-wave) and are snapped to zero.
        lateral_snap_mm = 0.01
        for axis in (0, 1):
            if abs(centre_mm[axis]) < lateral_snap_mm:
                centre_mm[axis] = 0.0
        if any(centre_mm[axis] != 0.0 for axis in (0, 1)):
            print(
                "  Warning: Significant lateral phase-centre offset; replicated port sheets may no "
                "longer lie exactly on the flush board edge."
            )
        centre = np.array(centre_mm)
        print(f"  Applied phase centre: {centre_mm} mm")
        print(f"  Plane-wise Z: phi=0 -> {z_phi0:.3f} mm, phi=90 -> {z_phi90:.3f} mm (astigmatism)")
        print(
            f"  Weighted RMS phase error over the cap: {rms_phase_error_deg(np.zeros(3)):.2f} deg "
            f"(origin) -> {rms_phase_error_deg(centre):.2f} deg (fitted centre)"
        )

        source_modeler = app.modeler
        _delete_coordinate_system(source_modeler, cs_name)
        print(f"  Saving '{cs_name}' coordinate system to unit cell...")
        source_modeler.create_coordinate_system(
            origin=[f"{c}mm" for c in centre_mm],
            reference_cs="Global",
            name=cs_name,
        )
        source_modeler.set_working_coordinate_system("Global")
        app.save_project()
        print("  Unit cell design saved.")
        return centre_mm

    def _unit_cell_far_field(self, sphere_name, expressions, freq_ghz, app=None):
        """
        Reads far-field expressions on `sphere_name` from the unit-cell solution: first from the
        setup's frequency sweeps at `freq_ghz`, then from LastAdaptive at the nearest adaptive
        frequency. Solves the setup once if no far-field data is available.
        """
        app = app or self.source_design_app
        setup_names = list(app.setup_names)
        if not setup_names:
            raise RuntimeError("Unit-cell design has no analysis setup to take the far field from.")
        setup_name = getattr(app, "active_setup", None) or setup_names[0]
        candidates = [
            (f"{setup_name} : {sweep}", [f"{freq_ghz}GHz"])
            for sweep in app.get_setup(setup_name).get_sweep_names()
        ] + [(f"{setup_name} : LastAdaptive", ["All"])]

        for attempt in range(2):
            for solution, freqs in candidates:
                try:
                    solution_data = app.post.get_solution_data(
                        expressions=expressions,
                        setup_sweep_name=solution,
                        variations={"Theta": ["All"], "Phi": ["All"], "Freq": freqs},
                        primary_sweep_variable="Theta",
                        report_category="Far Fields",
                        context=sphere_name,
                    )
                except Exception:
                    solution_data = None
                if not solution_data:
                    continue
                columns = _far_field_columns(solution_data, expressions)
                freq_column = columns["Freq"]
                if np.max(freq_column) > 1e6:  # reported in Hz
                    freq_column = freq_column / 1e9
                nearest = freq_column[np.argmin(np.abs(freq_column - freq_ghz))]
                keep = np.isclose(freq_column, nearest)
                print(f"  Far-field data: '{solution}' at {nearest:.3f} GHz")
                return {name: values[keep] for name, values in columns.items()}
            if attempt == 0:
                print(f"  No far-field data available; solving unit-cell setup '{setup_name}'...")
                app.analyze_setup(setup_name)
        raise RuntimeError("Could not obtain unit-cell far-field data for the phase-centre fit.")

    def close(self):
        """Releases the connection to the live AEDT session, which keeps running."""
        if self.desktop_session:
            self.desktop_session.release_desktop(close_projects=False, close_on_exit=False)
            print("AEDT desktop session connection released.")
            self.desktop_session = None


def load_layout_file(layout_path):
    """Loads a layout ('metadata' and 'elements') from a YAML (.yaml, .yml) or JSON file."""
    layout_path = os.path.abspath(layout_path)
    with open(layout_path, "r", encoding="utf-8") as f:
        if os.path.splitext(layout_path)[1].lower() == ".json":
            return json.load(f)
        import yaml

        return yaml.safe_load(f)


def main(argv=None):
    """Command-line entry point (installed as `hfss_array_builder`)."""
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Builds the MIMO array of a layout file in HFSS, attached to the running AEDT GUI "
            "session. The layout names its unit cells (metadata.primary_unit_cell, unit_cells) "
            "and modelling options."
        )
    )
    parser.add_argument("project_path", help="Path to the AEDT project (.aedt) with the unit cells")
    parser.add_argument("layout_path", help="Layout YAML or JSON file (layouts/generate_*.py)")
    parser.add_argument(
        "--output-design",
        default=None,
        metavar="NAME",
        help="Name of the array design to create (default: <primary_unit_cell><LayoutFileStem>)",
    )
    parser.add_argument(
        "--use-existing-cs",
        "--use-existing-phase-centre",
        "--use-existing-phase-center",
        action="store_true",
        default=None,
        help="Use the unit cells' existing PhaseCentreCS without asking (instead of a fit)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", default=None, help="Replace an existing array design"
    )
    simulate = parser.add_mutually_exclusive_group()
    simulate.add_argument(
        "--simulate",
        dest="simulate",
        action="store_true",
        default=None,
        help="Solve without asking",
    )
    simulate.add_argument(
        "--no-simulate", dest="simulate", action="store_false", help="Do not solve, do not ask"
    )
    args = parser.parse_args(argv)

    try:
        builder = MimoHfssBuilder(args.project_path, args.layout_path, args.output_design)
        print(
            f"Layout {args.layout_path}: primary unit cell '{builder.source_design_name}', "
            f"array design '{builder.target_design_name}'"
        )
        builder.connect_desktop()
        builder.synthesize_array_in_hfss(
            use_existing_cs=args.use_existing_cs,
            overwrite=args.overwrite,
            simulate=args.simulate,
        )
    except BuildAbortedError as e:
        print(f"[INFO] {e}")
        return 0
    except LayoutError as e:
        print(f"[ERROR] {e}")
        return 1
    except (NativeAedtRefusedError, NoLiveSessionError) as e:
        print(f"[ERROR] {e}")
        return 2
    finally:
        if "builder" in locals():
            builder.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
