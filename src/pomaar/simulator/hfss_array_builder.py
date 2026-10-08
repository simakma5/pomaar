#!/usr/bin/env python3
"""
Module 1: HFSS Full-Wave Array Synthesis.
Programmatically builds a planar MIMO antenna array layout on a single PCB substrate in HFSS
using a template-based, boolean-driven assembly workflow.
"""

import math
import os
import re
import sys

from ansys.aedt.core import Desktop, Hfss
from ansys.aedt.core.generic.constants import Axis, Gravity
from ansys.aedt.core.generic.data_handlers import _dict2arg
from ansys.aedt.core.generic.general_methods import grpc_active_sessions
import numpy as np

SPEED_OF_LIGHT_MM_S = 2.99792458e11
DEFAULT_CENTRE_FREQUENCY_GHZ = 79.0

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
# power normalisation); directivity is used so that XPD does not depend on port mismatch settings.
# With 'Linear' polarization (slant 0) and an unrotated CS, CoPolar/CrossPolar are Ludwig-3 X/Y
# referenced to the sphere CS x-axis: XPD of an x-polarised (H) port is Co - Cross, and of a
# y-polarised (V) port Cross - Co.
XPD_X_POL_EXPRESSION = "dB(DirCoPolar) - dB(DirCrossPolar)"
XPD_Y_POL_EXPRESSION = "dB(DirCrossPolar) - dB(DirCoPolar)"


class BuildAbortedError(Exception):
    """Raised when the user declines to continue at an interactive prompt."""


class NativeAedtRefusedError(RuntimeError):
    """Raised when connecting from a Linux host outside the AEDT container."""


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


def _parse_frequency_ghz(value):
    """Parses an AEDT frequency string (e.g. '77GHz', '500MHz') to GHz. Bare numbers are GHz."""
    if value is None or value == "":
        return None
    match = re.match(r"\s*([-+\d.eE]+)\s*([a-zA-Z]*)", str(value))
    if not match:
        return None
    number, unit = float(match.group(1)), match.group(2).lower()
    scale = {"": 1.0, "ghz": 1.0, "mhz": 1e-3, "khz": 1e-6, "hz": 1e-9, "thz": 1e3}
    return number * scale.get(unit, 1.0)


def _setup_frequency_ghz(app):
    """Returns the representative (centre) frequency of the first analysis setup in GHz, or None."""
    setup_names = list(app.setup_names)
    if not setup_names:
        return None
    setup = app.get_setup(getattr(app, "active_setup", None) or setup_names[0])
    props = setup.props

    single = _parse_frequency_ghz(props.get("Frequency"))
    if single:
        return single

    adaptive = props.get("MultipleAdaptiveFreqsSetup", {})
    if isinstance(adaptive, dict):
        if "Low" in adaptive and "High" in adaptive:
            low = _parse_frequency_ghz(adaptive["Low"])
            high = _parse_frequency_ghz(adaptive["High"])
            if low and high:
                return (low + high) / 2.0
        if isinstance(adaptive.get("AdaptAt"), list):
            freqs = [_parse_frequency_ghz(f.get("Frequency")) for f in adaptive["AdaptAt"]]
        else:
            freqs = [_parse_frequency_ghz(key) for key in adaptive.keys()]
        freqs = [f for f in freqs if f]
        if freqs:
            return sum(freqs) / len(freqs)

    for sweep_name in setup.get_sweep_names():
        sweep_props = getattr(setup.get_sweep(sweep_name), "props", {})
        start = _parse_frequency_ghz(sweep_props.get("RangeStart"))
        stop = _parse_frequency_ghz(sweep_props.get("RangeEnd"))
        if start and stop:
            return (start + stop) / 2.0
    return None


def _length_mm(app, expression):
    """
    Converts a length expression to mm. Plain literals are parsed directly so that PyAEDT does not
    create its temporary 'pyaedt_evaluator' variable in the design.
    """
    match = re.fullmatch(
        r"\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*(mm|um|nm|cm|m|mil|in)?\s*", str(expression)
    )
    if match:
        scale = {None: 1.0, "mm": 1.0, "um": 1e-3, "nm": 1e-6, "cm": 10.0, "m": 1e3,
                 "mil": 0.0254, "in": 25.4}
        return float(match.group(1)) * scale[match.group(2)]
    return float(app.evaluate_expression(expression)) * 1000.0


def _port_name_from_sheet(sheet_name):
    """Maps a replicated port sheet name (e.g. 'PortSheet_Tx_1_H') to its port ('Port_Tx_1_H')."""
    suffix = sheet_name.replace("PortSheet", "")
    if suffix.startswith("_"):
        suffix = suffix[1:]
    return f"Port_{suffix}"


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

    def __init__(
        self,
        project_path,
        source_design_name,
        target_design_name=None,
        pcb_margin_mm=0.0,
        grpc_port=None,
        non_graphical=True,
        centre_frequency_ghz=None,
        bandwidth_ghz=4.0,
        phase_reference_mode="single",
    ):
        self.project_path = os.path.abspath(project_path)
        self.source_design_name = source_design_name
        self.target_design_name = (
            target_design_name if target_design_name else f"{source_design_name}MimoArray"
        )
        self.pcb_margin_mm = pcb_margin_mm
        self.grpc_port = grpc_port
        self.non_graphical = non_graphical
        self.centre_frequency_ghz = centre_frequency_ghz
        self.bandwidth_ghz = bandwidth_ghz
        self.phase_reference_mode = phase_reference_mode
        self.port_to_cs_mapping = {}
        self.element_phase_centres = []

        self.desktop_session = None
        self.source_design_app = None
        self.target_design_app = None
        self.is_new_desktop = False
        self.is_new_design = True

    def connect_desktop(self):
        """Starts or connects to an existing AEDT session via gRPC."""
        # On a Linux host, PyAEDT cannot see a containerised AEDT's gRPC socket (separate network
        # namespace) and silently launches a native AEDT instead, which crashes there. Run the
        # builder inside the container (`ansys-vnc exec ...`) or set POMAAR_ALLOW_NATIVE_AEDT=1.
        if (
            sys.platform.startswith("linux")
            and not os.path.exists("/run/.containerenv")
            and os.environ.get("POMAAR_ALLOW_NATIVE_AEDT") != "1"
        ):
            raise NativeAedtRefusedError(
                "Refusing to connect from outside the AEDT container: run this via "
                "`ansys-vnc exec`, or set POMAAR_ALLOW_NATIVE_AEDT=1 for a native AEDT install."
            )
        if self.grpc_port is None:
            # Attach to the open GUI session, whose gRPC port is the first free one from 50051
            graphical_ports = grpc_active_sessions(non_graphical=False)
            self.grpc_port = graphical_ports[0] if graphical_ports else 50051
        print(f"Connecting to or starting AEDT on port {self.grpc_port}...")
        try:
            self.desktop_session = Desktop(
                port=self.grpc_port,
                new_desktop=False,
                non_graphical=self.non_graphical,
                close_on_exit=False,
            )
            print("Connected to existing active AEDT session.")
            self.is_new_desktop = False
        except Exception as err:
            print(f"Failed to connect to active session ({err}). Starting new session...")
            self.desktop_session = Desktop(
                port=self.grpc_port,
                new_desktop=True,
                non_graphical=self.non_graphical,
                close_on_exit=True,
            )
            self.is_new_desktop = True

    def calculate_default_coplanar_layout(
        self,
        transmitter_count=16,
        receiver_count=16,
        operating_frequency_ghz=79.0,
        subarray_spacing_mm=10.0,
    ):
        """
        Calculates coordinates for a coplanar MIMO array on the z=0 plane.

        The Rx array is a dense ULA along the X-axis (centred at y=0, z=0) with 0.5 lambda spacing.
        The Tx array is a sparse ULA along the X-axis (centred at y=subarray_spacing, z=0)
        with spacing = receiver_count * receiver_spacing.
        """
        # Speed of light in vacuum (mm/s)
        speed_of_light_mm_s = 2.99792458e11
        wavelength_mm = speed_of_light_mm_s / (operating_frequency_ghz * 1e9)
        receiver_spacing_mm = 0.5 * wavelength_mm
        transmitter_spacing_mm = receiver_count * receiver_spacing_mm

        elements = []

        # Rx elements along X-axis at y=0
        rx_offset_x = (receiver_count - 1) * receiver_spacing_mm / 2.0
        for rx_idx in range(receiver_count):
            x_pos = rx_idx * receiver_spacing_mm - rx_offset_x
            elements.append(
                {
                    "label": f"Rx_{rx_idx + 1}",
                    "pos": [x_pos, 0.0, 0.0],
                    "role": "Rx",
                    "polarization": "v",
                    "yaw": 0.0,
                }
            )

        # Tx elements along X-axis at y=subarray_spacing
        tx_offset_x = (transmitter_count - 1) * transmitter_spacing_mm / 2.0
        for tx_idx in range(transmitter_count):
            x_pos = tx_idx * transmitter_spacing_mm - tx_offset_x
            elements.append(
                {
                    "label": f"Tx_{tx_idx + 1}",
                    "pos": [x_pos, subarray_spacing_mm, 0.0],
                    "role": "Tx",
                    "polarization": "v",
                    "yaw": 180.0,
                }
            )

        print(f"Generated coplanar layout at {operating_frequency_ghz} GHz:")
        print(f"  Rx Elements: {receiver_count} (spacing={receiver_spacing_mm:.2f} mm)")
        print(f"  Tx Elements: {transmitter_count} (spacing={transmitter_spacing_mm:.2f} mm)")
        return elements

    def synthesize_array_in_hfss(
        self,
        elements,
        operating_frequency_ghz=None,
        setup_results=None,
        run_simulation=None,
        use_existing_cs=None,
        overwrite=None,
    ):
        """
        Synthesizes the MIMO array in HFSS using a naming-convention boolean assembly:
        1. Classifies objects in the unit cell (Global Layers, Port Sheets, Active Copper, Dummy Cutouts).
        2. Sizes and draws the global continuous substrate and ground layers.
        3. Copies the source templates into the target design exactly ONCE.
        4. Replicates active structures, port sheets, and dummy solids to all coordinates using local duplicate command.
        5. Applies boolean operations using the dummy solids (e.g. Subtract_L2_Ground).
        6. Assigns lumped port excitations to the replicated port sheets.
        7. Creates radiation airbox, sets up simulation sweep, asks user whether to set up results, and prompts to launch simulation ('Analyze All').

        Raises
        ------
        BuildAbortedError
            If the user declines to continue at one of the interactive prompts.
        """
        if isinstance(elements, dict):
            layout_metadata = elements.get("metadata", {})
            elements_list = elements.get("elements", [])
        else:
            layout_metadata = {}
            elements_list = elements

        if "bandwidth_ghz" in layout_metadata:
            self.bandwidth_ghz = float(layout_metadata["bandwidth_ghz"])

        if not self.desktop_session:
            self.connect_desktop()

        print(f"Loading project: {self.project_path}")
        self.source_design_app = Hfss(
            project=self.project_path,
            design=self.source_design_name,
            new_desktop=False,
            close_on_exit=False,
        )

        # Centre frequency precedence: explicit argument > constructor > layout metadata >
        # unit-cell analysis setup > default.
        if operating_frequency_ghz is not None:
            self.centre_frequency_ghz = float(operating_frequency_ghz)
        elif self.centre_frequency_ghz is None:
            if "center_frequency_ghz" in layout_metadata:
                self.centre_frequency_ghz = float(layout_metadata["center_frequency_ghz"])
            else:
                self.centre_frequency_ghz = (
                    _setup_frequency_ghz(self.source_design_app) or DEFAULT_CENTRE_FREQUENCY_GHZ
                )
        operating_frequency_ghz = self.centre_frequency_ghz
        print(f"Centre frequency: {operating_frequency_ghz:.3f} GHz")

        source_modeler = self.source_design_app.modeler
        all_object_names = list(source_modeler.object_names)

        # --- STEP 1: Classify Objects ---
        print("\nScanning and classifying objects in the unit-cell design:")
        global_layers = {}
        active_elements = []
        port_sheets = []
        dummy_solids = {}

        for name in all_object_names:
            obj = source_modeler.get_object_from_name(name)
            if not obj:
                continue

            # Check for boolean dummy solids first (e.g. Subtract_L2_Ground)
            if name.startswith("Subtract_") or name.startswith("Unite_"):
                tokens = name.split("_")
                operation = tokens[0]
                target_solid = "_".join(tokens[1:])
                key = (operation, target_solid)
                if key not in dummy_solids:
                    dummy_solids[key] = []
                dummy_solids[key].append(name)
                print(f"  [Dummy {operation}]   {name} -> target: {target_solid}")
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
                global_layers[name] = {
                    "material": obj.material_name,
                    "z_min": z_min,
                    "z_max": z_max,
                }
                tag = "[Global Layer]    " if is_substrate else "[Global Ground]   "
                print(
                    f"  {tag} {name} ({obj.material_name}, Z=[{z_min:.2f}, {z_max:.2f}] mm)"
                )

            elif "PortSheet" in name or name.startswith("Port_"):
                port_sheets.append(name)
                print(f"  [Port Sheet]       {name}")

            else:
                active_elements.append(name)
                print(f"  [Active Element]   {name}")

        run_fit = False
        offset = [0.0, 0.0, 0.0]

        cs_map = {cs.name: cs for cs in source_modeler.coordinate_systems}
        if "PhaseCentreCS" in cs_map:
            cs_props = cs_map["PhaseCentreCS"].props
            try:
                offset = [
                    round(_length_mm(self.source_design_app, cs_props[key]), 3)
                    for key in ("OriginX", "OriginY", "OriginZ")
                ]
            except Exception as e:
                print(
                    f"  Warning: Could not evaluate PhaseCentreCS origin ({e}); using zero offset."
                )
            print(f"\n[INFO] PhaseCentreCS found in unit-cell design. Origin offset: {offset} mm.")
            if not _confirm(
                "Do you want to use the existing PhaseCentreCS?", preset=use_existing_cs
            ):
                run_fit = True
        else:
            print("\n" + "=" * 80)
            print("[WARNING] PhaseCentreCS was NOT found in the unit-cell design!")
            print("=" * 80 + "\n")
            if use_existing_cs:
                print("[INFO] Proceeding without PhaseCentreCS offset [0.0, 0.0, 0.0] mm.")
            elif _confirm("Do you want to fit the phase centre from the unit-cell far field?"):
                run_fit = True
            elif _confirm("Proceed without PhaseCentreCS offset?"):
                print("[INFO] Proceeding with zero offset [0.0, 0.0, 0.0] mm.")
            else:
                raise BuildAbortedError("Aborted by user: no phase centre available.")

        if run_fit:
            offset = self.fit_unit_cell_phase_centre()

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
            ):
                raise BuildAbortedError("Aborted by user to prevent overwriting existing design.")

            print(f"Design '{self.target_design_name}' already exists. Reusing and clearing it...")
            self.target_design_app = Hfss(
                project=project_name,
                design=self.target_design_name,
                new_desktop=False,
                close_on_exit=False,
            )
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
                self.target_design_app = Hfss(
                    project=project_name,
                    design=self.target_design_name,
                    new_desktop=False,
                    close_on_exit=False,
                )

        # Reconnect source_design_app connection to make sure it points to the unit cell design
        self.source_design_app = Hfss(
            project=project_name,
            design=self.source_design_name,
            new_desktop=False,
            close_on_exit=False,
        )

        # Sync design variables from source design to target design
        # This prevents copied objects in the target design from retaining stale/cached design variable dimensions.
        try:
            print("Syncing design variables from source to target design...")
            for (
                var_name,
                var_expr,
            ) in self.source_design_app.variable_manager.design_variables.items():
                self.target_design_app[var_name] = var_expr
        except Exception as e:
            print(f"  Warning: Failed to sync design variables ({e})")

        target_modeler = self.target_design_app.modeler

        # Measure unit-cell bounding box relative to PhaseCentreCS (offset)
        # Note: We measure exclusively from active elements and port sheets so that
        # the measurement reflects the true element reach (L_feed) and is never corrupted
        # by pre-existing or oversized substrate boundaries.
        uc_x_extents = []
        uc_y_extents = []
        for obj_name in list(global_layers.keys()) + active_elements + port_sheets:
            try:
                obj = self.source_design_app.modeler.get_object_from_name(obj_name)
                if obj:
                    bb = obj.bounding_box
                    uc_x_extents.extend(
                        [abs(float(bb[0]) - offset[0]), abs(float(bb[3]) - offset[0])]
                    )
                    uc_y_extents.extend(
                        [abs(float(bb[1]) - offset[1]), abs(float(bb[4]) - offset[1])]
                    )
            except Exception:
                pass

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
        self.target_design_app["unitCellHalfWidth"] = f"{unit_cell_extent_x:.4f}mm"
        self.target_design_app["unitCellHalfLength"] = f"{unit_cell_extent_y:.4f}mm"
        self.target_design_app["airboxClearance"] = f"{airbox_clearance_mm:.4f}mm"
        self.target_design_app["pcbMargin"] = f"{self.pcb_margin_mm:.4f}mm"

        # Inject layout variables from layout metadata
        layout_vars = layout_metadata.get("variables", {})
        for var_name, var_expr in layout_vars.items():
            self.target_design_app[var_name] = str(var_expr)

        # --- STEP 3: Size & Draw Global Contiguous Board ---
        # Sort layers by prefix (L1, L12, L2...) to ensure they are created in stackup order
        sorted_layers = sorted(global_layers.keys(), key=lambda name: name.split("_")[0])
        ground_layer_names = []

        # Track overall board bounds for numerical fallback and Z calculation
        overall_x_min = None
        overall_x_max = None
        overall_y_min = None
        overall_y_max = None

        # Y-reach of the active structures (feedlines and ports) is layer-independent
        active_y_extents = []
        for active_name in active_elements + port_sheets:
            try:
                act_obj = self.source_design_app.modeler.get_object_from_name(active_name)
                act_bbox = act_obj.bounding_box
                active_y_extents.extend([float(act_bbox[1]), float(act_bbox[4])])
            except Exception:
                pass

        for layer_name in sorted_layers:
            layer_obj = self.source_design_app.modeler.get_object_from_name(layer_name)
            bbox = layer_obj.bounding_box
            x_min_uc, x_max_uc = float(bbox[0]), float(bbox[3])

            y_min_active_uc = min(active_y_extents) if active_y_extents else None
            y_max_active_uc = max(active_y_extents) if active_y_extents else None

            if y_min_active_uc is None:
                y_min_active_uc = float(bbox[1])
            if y_max_active_uc is None:
                y_max_active_uc = float(bbox[4])

            corners = [
                (x_min_uc, y_min_active_uc),
                (x_max_uc, y_min_active_uc),
                (x_min_uc, y_max_active_uc),
                (x_max_uc, y_max_active_uc),
            ]

            all_x_glob = []
            all_y_glob = []
            for element in elements_list:
                pos = element.get("position", element.get("pos", [0.0, 0.0, 0.0]))
                yaw_deg = float(
                    element.get("yaw", element.get("rotation_yaw", element.get("rotation", 0.0)))
                )
                rad = math.radians(yaw_deg)
                cos_val = math.cos(rad)
                sin_val = math.sin(rad)

                for cx, cy in corners:
                    rx = cx * cos_val - cy * sin_val
                    ry = cx * sin_val + cy * cos_val
                    all_x_glob.append(rx + pos[0])
                    all_y_glob.append(ry + pos[1])

            global_x_min = min(all_x_glob) - self.pcb_margin_mm
            global_x_max = max(all_x_glob) + self.pcb_margin_mm
            global_y_min = min(all_y_glob) - self.pcb_margin_mm
            global_y_max = max(all_y_glob) + self.pcb_margin_mm

            if overall_x_min is None or global_x_min < overall_x_min:
                overall_x_min = global_x_min
            if overall_x_max is None or global_x_max > overall_x_max:
                overall_x_max = global_x_max
            if overall_y_min is None or global_y_min < overall_y_min:
                overall_y_min = global_y_min
            if overall_y_max is None or global_y_max > overall_y_max:
                overall_y_max = global_y_max

        # Determine arrayBoardWidth and arrayBoardLength formulas or values
        # We use arrayBoardWidth / arrayBoardLength for the synthesized board layers so that
        # any internal unit-cell variables (like boardWidth / boardLength / sub_L) retain their
        # single-element values without double-expanding the feedlines.
        board_config = layout_metadata.get("board", {})
        width_formula = board_config.get("width_formula")
        length_formula = board_config.get("length_formula")

        board_width_calc = overall_x_max - overall_x_min if overall_x_max is not None else 10.0
        board_height_calc = overall_y_max - overall_y_min if overall_y_max is not None else 10.0

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
                        f"{z_min - offset[2]:.4f}mm",
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
                        f"  Creating ground plane sheet '{layer_name}' at Z={z_min - offset[2]:.2f} mm"
                    )
                    target_modeler.create_rectangle(
                        orientation="XY",
                        origin=[
                            "-arrayBoardWidth / 2",
                            "-arrayBoardLength / 2",
                            f"{z_min - offset[2]:.4f}mm",
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
                            f"{z_min - offset[2]:.4f}mm",
                        ],
                        sizes=["arrayBoardWidth", "arrayBoardLength", f"{thickness:.4f}mm"],
                        name=layer_name,
                        material=material,
                    )

        # --- STEP 4: Replicate Elements via Dedicated Coordinate Systems ---
        all_replicate_sources = active_elements + port_sheets
        for dummy_list in dummy_solids.values():
            all_replicate_sources.extend(dummy_list)

        replicated_dummy_mapping = {key: [] for key in dummy_solids.keys()}
        replicated_ports_list = []
        element_phase_centres = []
        element_phase_centre_exprs = []
        element_placements = []  # (label, [x, y, z] mm, yaw deg) for the post-build check
        port_to_cs_mapping = {}

        # Optimization: Copy the templates from the source design to the target design exactly ONCE.
        print("\nCopying template geometries to target design...")
        target_modeler.set_working_coordinate_system("Global")
        self.target_design_app.copy_solid_bodies_from(
            design=self.source_design_app,
            assignment=all_replicate_sources,
            no_vacuum=False,
            no_pec=False,
            include_sheets=True,
        )

        # Clear any copied boundaries/excitations to strip wave ports from template geometries
        try:
            self.target_design_app.oboundary.DeleteAllBoundaries()
            self.target_design_app.oboundary.DeleteAllExcitations()
            print("  Stripped lingering wave port assignments from template sheets.")
        except Exception as e:
            print(f"  Warning: Failed to clear boundaries/excitations ({e})")

        self.target_design_app.set_active_design(self.target_design_name)

        print("Replicating structures to array grid with element coordinate systems...")
        for element in elements_list:
            raw_label = element.get("label", element.get("name", "Element"))
            pos = element.get("position", element.get("pos", [0.0, 0.0, 0.0]))
            pos_expr = element.get("position_expression", element.get("pos_expr", None))
            element_yaw = float(
                element.get("yaw", element.get("rotation_yaw", element.get("rotation", 0.0)))
            )
            pol = str(element.get("polarization", element.get("pol", ""))).strip().upper()

            # Append polarization suffix if specified and not already present in the label
            if pol and not raw_label.upper().endswith(f"_{pol}"):
                label = f"{raw_label}_{pol}"
            else:
                label = raw_label

            cs_name = f"CS_{label}"
            print(f"\n  Setting up element {label} (yaw={element_yaw:.1f} deg)...")

            # Track the element phase centre in array coordinates. The whole stack is shifted by
            # -offset[2] in Z, which places the unit-cell phase centre exactly at pos[2].
            z_pos_val = float(pos[2]) if len(pos) > 2 else 0.0
            element_phase_centres.append([float(pos[0]), float(pos[1]), z_pos_val])

            # Element CS sits at the element phase centre (parametric where available). Its origin Z
            # does not affect placement: only its axes and vertical rotation axis are used below.
            if pos_expr:
                cs_origin = [str(pos_expr[0]), str(pos_expr[1]), f"{z_pos_val:.4f}mm"]
                element_phase_centre_exprs.append([str(pos_expr[0]), str(pos_expr[1])])
            else:
                cs_origin = [f"{pos[0]:.4f}mm", f"{pos[1]:.4f}mm", f"{z_pos_val:.4f}mm"]
                element_phase_centre_exprs.append([f"{pos[0]:.4f}mm", f"{pos[1]:.4f}mm"])

            # Calculate pointing vectors for the CS rotation around Z
            rad = math.radians(element_yaw)
            cos_val = math.cos(rad)
            sin_val = math.sin(rad)
            x_pointing = [cos_val, sin_val, 0.0]
            y_pointing = [-sin_val, cos_val, 0.0]

            # Recreate coordinate system if it exists
            _delete_coordinate_system(target_modeler, cs_name)

            print(f"  Creating relative coordinate system '{cs_name}' at origin {cs_origin}...")
            target_modeler.create_coordinate_system(
                origin=cs_origin,
                reference_cs="Global",
                name=cs_name,
                mode="axis",
                x_pointing=x_pointing,
                y_pointing=y_pointing,
            )
            # Creating a CS makes it the working CS. All placement operations below are defined
            # in Global coordinates (rotation about the global Z axis, then a global translation),
            # so the working CS must be reset; otherwise the rotation pivots about the element CS
            # and the move vector is read in its rotated axes.
            target_modeler.set_working_coordinate_system("Global")
            element_placements.append(
                (label, [float(pos[0]), float(pos[1]), z_pos_val], element_yaw)
            )

            # Duplicate all templates locally using duplicate_along_line with a dummy Z offset of 1.0 mm.
            success, pasted_names = target_modeler.duplicate_along_line(
                assignment=all_replicate_sources,
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

            # Rename duplicated objects, assign to element CS, and track ports/dummy solids
            renamed_objs = []
            for pasted_name in pasted_names:
                # Strip numerical suffixes added by AEDT duplicate if any (e.g. L1_Patch_1 -> L1_Patch)
                base_name = pasted_name
                for source_name in all_replicate_sources:
                    if pasted_name.startswith(source_name):
                        base_name = source_name
                        break

                new_name = f"{base_name}_{label}"
                obj = target_modeler.get_object_from_name(pasted_name)
                obj.name = new_name
                renamed_objs.append(new_name)

                # Track port sheets
                if base_name in port_sheets:
                    replicated_ports_list.append(new_name)
                    port_to_cs_mapping[_port_name_from_sheet(new_name)] = cs_name

                # Track dummy solids
                for key, dummy_src_list in dummy_solids.items():
                    if base_name in dummy_src_list:
                        replicated_dummy_mapping[key].append(new_name)

            # Rotate element if needed (about the global Z axis through the unit-cell origin)
            if element_yaw != 0.0:
                target_modeler.rotate(
                    assignment=renamed_objs,
                    axis=Axis.Z,
                    angle=element_yaw,
                )

            # Translate so that the rotated phase centre R*offset lands on the layout position:
            # p -> R(p) + pos - R(offset) = R(p - offset) + pos
            rot_dx = offset[0] * cos_val - offset[1] * sin_val
            rot_dy = offset[0] * sin_val + offset[1] * cos_val
            rot_dz = offset[2]

            # In Z, the element is translated by (pos[2] - offset[2]), matching the -offset[2]
            # shift of the board layers so the stack-up is preserved exactly.
            z_shift_str = f"{z_pos_val - rot_dz:.4f}mm"

            # Move element to final position
            if pos_expr:
                # Parenthesised so that negative offsets and unary minus in the layout
                # expression never produce ambiguous terms such as "-r - -0.005mm"
                x_move = (
                    f"({pos_expr[0]}) - ({rot_dx:.4f}mm)"
                    if abs(rot_dx) > 1e-4
                    else str(pos_expr[0])
                )
                y_move = (
                    f"({pos_expr[1]}) - ({rot_dy:.4f}mm)"
                    if abs(rot_dy) > 1e-4
                    else str(pos_expr[1])
                )
                move_vector = [x_move, y_move, z_shift_str]
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

        # Reset working coordinate system to Global
        target_modeler.set_working_coordinate_system("Global")
        self.element_phase_centres = element_phase_centres
        self.port_to_cs_mapping = port_to_cs_mapping

        # Clean up the original template geometries at the origin of the target design
        print("\nCleaning up template geometries...")
        target_modeler.delete(all_replicate_sources)

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
        for port_sheet in replicated_ports_list:
            port_name = _port_name_from_sheet(port_sheet)

            print(f"  Assigning wave port to sheet '{port_sheet}' referencing '{default_ground}'")
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
        port_names = sorted(_port_name_from_sheet(sheet) for sheet in replicated_ports_list)
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
        # PhaseCentreCS sits at the centroid of the element phase centres with global axes. It is
        # the CS of all far-field spheres (z = boresight, x = co-polar reference) and the single
        # far-field phase reference. Lateral coordinates stay parametric when the layout is.
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

        # Set Far-field Phase Reference for SBR+ Antenna Linking
        self.set_far_field_phase_reference(mode=self.phase_reference_mode)

        # --- STEP 7: Create Radiation Airbox ---
        print("\nCreating radiation boundary Airbox...")
        # Get overall board bounds in Z direction
        all_z_coords = []
        for layer_info in global_layers.values():
            all_z_coords.extend([layer_info["z_min"], layer_info["z_max"]])
        overall_z_min = (min(all_z_coords) if all_z_coords else 0.0) - offset[2]
        overall_z_max = (max(all_z_coords) if all_z_coords else 1.0) - offset[2]
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
            element_placements, active_elements + port_sheets, offset, replicated_ports_list
        )
        if placement_errors:
            self.target_design_app.save_project()
            raise RuntimeError(
                "Post-build placement check failed (project saved for inspection):\n  "
                + "\n  ".join(placement_errors)
            )

        # Configure simulation setup and sweep by inheriting from the element design
        self.configure_simulation_setup()

        # As the last step of building, ask the user whether to set up results (far-field sphere & post-processing reports)
        if _confirm(
            "\nDo you want to set up results (far-field spheres & post-processing reports)?",
            preset=setup_results,
        ):
            self.create_post_processing_reports()
        else:
            print("\n[INFO] Skipping results setup as requested.")

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

        # As the very last step, ask the user if they want to launch the simulation ('Analyze All')
        if _confirm(
            "\nDo you want to launch the simulation ('Analyze All')?", preset=run_simulation
        ):
            print("\nLaunching HFSS simulation ('Analyze All')...")
            self.target_design_app.analyze()
        else:
            print("\n[INFO] Skipping simulation execution. Model synthesis is complete.")

    def verify_placement(
        self, element_placements, template_names, offset, port_sheet_names, tolerance_mm=1e-3
    ):
        """
        Checks the built geometry against the layout, independently of how AEDT applied the
        placement operations.

        1. Every replicated template object must satisfy p -> R(yaw) (p - offset_xy) + pos_xy and
           z -> z - offset_z + pos_z, compared via the vertex centroid (exact under any rotation).
        2. Every wave port sheet must lie on a lateral face of the Airbox (exterior boundary).

        Returns a list of human-readable violations (empty when everything is in place).
        """
        source_modeler = self.source_design_app.modeler
        target_modeler = self.target_design_app.modeler
        print("\nVerifying element placement and port positions...")

        def vertex_centroid(obj):
            positions = [vertex.position for vertex in obj.vertices]
            return np.mean(np.array(positions, dtype=float), axis=0) if positions else None

        errors = []
        templates = {}
        for name in template_names:
            centroid = vertex_centroid(source_modeler[name])
            if centroid is None:
                print(f"  Note: '{name}' has no vertices; its placement is not checked.")
            else:
                templates[name] = centroid

        checked = 0
        for label, pos, yaw in element_placements:
            cos_val, sin_val = math.cos(math.radians(yaw)), math.sin(math.radians(yaw))
            for name, template_centroid in templates.items():
                obj = target_modeler.get_object_from_name(f"{name}_{label}")
                if not obj:
                    errors.append(f"Missing replicated object '{name}_{label}'")
                    continue
                dx, dy = template_centroid[0] - offset[0], template_centroid[1] - offset[1]
                expected = np.array(
                    [
                        cos_val * dx - sin_val * dy + pos[0],
                        sin_val * dx + cos_val * dy + pos[1],
                        template_centroid[2] - offset[2] + pos[2],
                    ]
                )
                actual = vertex_centroid(obj)
                deviation = np.linalg.norm(actual - expected)
                checked += 1
                if deviation > tolerance_mm:
                    errors.append(
                        f"'{name}_{label}' is {deviation:.4f} mm off its layout position "
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
        Configures the adaptive mesh setup and frequency sweep in the target design,
        inheriting the configuration (single/multi-frequency/broadband meshing,
        frequency specifications, lambda refinement, passes, and all frequency sweeps)
        from the single-element design.

        First attempts to natively copy-paste the setup from the source design module.
        If native copy-paste is unavailable, recreates the setup and sweeps manually
        from the source design's properties. If no setup exists in the source design,
        creates a default broadband setup.
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

        target_setup_name = "ArraySetup"
        source_setups = list(self.source_design_app.setup_names) if self.source_design_app else []
        copied_successfully = False

        if source_setups:
            src_setup_name = getattr(self.source_design_app, "active_setup", None) or source_setups[0]
            print(
                f"\nInheriting analysis setup and frequency sweep from element design ('{src_setup_name}')..."
            )
            # Primary approach: Native CopyDrivenSetup / PasteDrivenSetup via oanalysis
            try:
                self.source_design_app.oanalysis.CopyDrivenSetup(src_setup_name)
                pasted_name = self.target_design_app.oanalysis.PasteDrivenSetup()
                if pasted_name:
                    if pasted_name != target_setup_name:
                        self.target_design_app.oanalysis.RenameSetup(pasted_name, target_setup_name)
                    copied_successfully = True
                    print(f"  Successfully copy-pasted and renamed setup to '{target_setup_name}'.")
            except Exception as err:
                print(
                    f"  Warning: Native CopyDrivenSetup/PasteDrivenSetup encountered an issue ({err}). Falling back to manual recreation..."
                )

            # Secondary approach: Manual recreation from source properties
            if not copied_successfully:
                try:
                    print(f"  Recreating setup '{target_setup_name}' from source properties...")
                    src_setup = self.source_design_app.get_setup(src_setup_name)
                    solve_type = str(src_setup.props.get("SolveType", "BroadBand")).lower()

                    setup = self.target_design_app.create_setup(name=target_setup_name)
                    setup.auto_update = False

                    max_passes = src_setup.props.get("MaximumPasses", 21)
                    max_delta_s = src_setup.props.get("MaxDeltaS", 0.02)

                    if solve_type == "broadband":
                        adaptive_freqs = src_setup.props.get("MultipleAdaptiveFreqsSetup", {})
                        low_f = adaptive_freqs.get("Low", f"{self.centre_frequency_ghz - self.bandwidth_ghz / 2.0}GHz")
                        high_f = adaptive_freqs.get("High", f"{self.centre_frequency_ghz + self.bandwidth_ghz / 2.0}GHz")
                        setup.enable_adaptive_setup_broadband(
                            low_frequency=low_f,
                            high_frequency=high_f,
                            max_passes=max_passes,
                            max_delta_s=max_delta_s,
                        )
                        print(f"  Configured Broadband adaptive mesh: Low={low_f}, High={high_f}, MaxPasses={max_passes}, MaxDeltaS={max_delta_s}")
                    elif solve_type == "multifrequency":
                        adaptive_freqs = src_setup.props.get("MultipleAdaptiveFreqsSetup", {})
                        freq_list = list(adaptive_freqs.keys()) if isinstance(adaptive_freqs, dict) else [f"{self.centre_frequency_ghz}GHz"]
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

                    # Recreate frequency sweeps
                    try:
                        src_sweeps = src_setup.get_sweep_names()
                    except Exception:
                        src_sweeps = []

                    for sw_name in src_sweeps:
                        sw_src = src_setup.get_sweep(sw_name)
                        if hasattr(sw_src, "props"):
                            sw_props = sw_src.props
                            range_type = sw_props.get("RangeType", "LinearCount")
                            sw_type = sw_props.get("Type", "Interpolating")
                            start_f = sw_props.get("RangeStart", f"{self.centre_frequency_ghz - self.bandwidth_ghz / 2.0}GHz")
                            end_f = sw_props.get("RangeEnd", f"{self.centre_frequency_ghz + self.bandwidth_ghz / 2.0}GHz")
                            save_fields = sw_props.get("SaveFields", True)

                            if range_type == "LinearCount":
                                count = sw_props.get("RangeCount", 401)
                                self.target_design_app.create_linear_count_sweep(
                                    setup=target_setup_name,
                                    unit="GHz" if ("GHz" in str(start_f) or "GHz" in str(end_f)) else "",
                                    start_frequency=start_f,
                                    stop_frequency=end_f,
                                    num_of_freq_points=count,
                                    name=sw_name,
                                    save_fields=save_fields,
                                    sweep_type=sw_type,
                                )
                                print(f"  Created {sw_type} sweep '{sw_name}': {start_f} to {end_f} ({count} points, SaveFields={save_fields})")
                            elif range_type == "LinearStep":
                                step = sw_props.get("RangeStep", "10MHz")
                                self.target_design_app.create_linear_step_sweep(
                                    setup=target_setup_name,
                                    unit="GHz" if ("GHz" in str(start_f) or "GHz" in str(end_f)) else "",
                                    start_frequency=start_f,
                                    stop_frequency=end_f,
                                    step_size=step,
                                    name=sw_name,
                                    save_fields=save_fields,
                                    sweep_type=sw_type,
                                )
                                print(f"  Created {sw_type} sweep '{sw_name}': {start_f} to {end_f} (step {step}, SaveFields={save_fields})")
                    copied_successfully = True
                except Exception as err2:
                    print(f"  Warning: Manual recreation failed ({err2}). Falling back to default setup...")

        # Fallback if source design has no setups or if recreation failed
        if not copied_successfully:
            start_freq = self.centre_frequency_ghz - self.bandwidth_ghz / 2.0
            end_freq = self.centre_frequency_ghz + self.bandwidth_ghz / 2.0
            print(
                f"\nConfiguring default broadband adaptive mesh setup '{target_setup_name}' ({start_freq:.2f} GHz - {end_freq:.2f} GHz)..."
            )
            setup = self.target_design_app.create_setup(name=target_setup_name)
            setup.enable_adaptive_setup_broadband(
                low_frequency=f"{start_freq:.4f}GHz",
                high_frequency=f"{end_freq:.4f}GHz",
                max_passes=21,
                max_delta_s=0.02,
            )
            sweep_name = "Sweep"
            print(
                f"Configuring interpolating frequency sweep '{sweep_name}' ({start_freq:.2f} GHz - {end_freq:.2f} GHz, 401 points)..."
            )
            self.target_design_app.create_linear_count_sweep(
                setup=target_setup_name,
                unit="GHz",
                start_frequency=start_freq,
                stop_frequency=end_freq,
                num_of_freq_points=401,
                name=sweep_name,
                save_fields=True,
                sweep_type="Interpolating",
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
            # Strip mode suffixes (e.g. "Port_Rx_1_V:1" -> "Port_Rx_1_V") and deduplicate
            port_names = sorted(list(set(p.split(":")[0] for p in raw_excitations if p)))
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
            s_name = "ArraySetup" if "ArraySetup" in target_setups else target_setups[0]
            try:
                s_obj = self.target_design_app.get_setup(s_name)
                sw_names = s_obj.get_sweep_names()
                sw_name = sw_names[0] if sw_names else "Sweep"
            except Exception:
                sw_name = "Sweep"
            setup_sweep = f"{s_name} : {sw_name}"
        else:
            setup_sweep = "ArraySetup : Sweep"

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

    def run_solve(self):
        """Triggers the HFSS simulation setup and solve using 'Analyze All'."""
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

        # Note: the setup is deliberately not reconfigured here, since that deletes existing
        # setups together with their solutions.
        print("Solving full-wave HFSS array design using 'Analyze All'...")
        self.target_design_app.analyze()

    def export_coupling_s_parameters(self, output_touchstone_path, setup_name="ArraySetup"):
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
            - "single" (default): Assigns a single coordinate system (array centroid 'PhaseCentreCS')
              to all excitations. This preserves the inter-element spatial phase differences (array manifold)
              necessary for MIMO radar SBR+ simulations.
            - "per_port" or "per-port": Assigns individual coordinate systems to each port
              (e.g. 'CS_<label>' for each element).
        cs_name : str, optional
            The coordinate system name to use when mode="single". Default is "PhaseCentreCS".
        """
        if not self.target_design_app:
            raise RuntimeError("Target array design not connected or synthesized yet.")

        if mode is None:
            mode = getattr(self, "phase_reference_mode", "single")
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
                print(f"  Warning: Failed to set single phase center to '{target_cs}' ({e})")
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
                print(f"  Warning: Failed to set per-port phase center ({e})")
        else:
            print(f"  Warning: Unknown phase reference mode '{mode}'. Skipping.")

    def fit_unit_cell_phase_centre(self, cap_deg=40.0, cs_name="PhaseCentreCS"):
        """
        Fits the unit-cell phase centre from its solved far field and stores it as `cs_name` in the
        unit-cell design.

        The co-polar Ludwig-3 component (L3X or L3Y, whichever dominates at boresight) is sampled
        over the cap theta <= cap_deg (all phi) at the centre frequency. A phase centre at d makes
        the far-field phase vary as arg E_co(u) = psi_0 + k u.d (u = unit direction), so d is the
        |E_co|-weighted least-squares solution of that linear model. Unlike a single-cut phase
        flatness optimisation, this constrains X, Y and Z and needs no extra solves. The separate
        phi = 0 / 90 plane fits are reported as a measure of astigmatism.

        Returns
        -------
        list of float
            Phase centre [x, y, z] in mm, in unit-cell global coordinates.
        """
        app = self.source_design_app
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
            data = self._unit_cell_far_field(sphere_name, ["rEL3X", "rEL3Y"], freq_ghz)
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

    def _unit_cell_far_field(self, sphere_name, expressions, freq_ghz):
        """
        Reads far-field expressions on `sphere_name` from the unit-cell solution: first from the
        setup's frequency sweeps at `freq_ghz`, then from LastAdaptive at the nearest adaptive
        frequency. Solves the setup once if no far-field data is available.
        """
        app = self.source_design_app
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
        """Safely releases the AEDT connection."""
        if self.desktop_session:
            if self.is_new_desktop:
                self.desktop_session.close_desktop()
                print("AEDT desktop session closed.")
            else:
                self.desktop_session.release_desktop(close_projects=False, close_on_exit=False)
                print("AEDT desktop session connection released.")
            self.desktop_session = None


def load_layout_file(layout_path):
    """
    Loads an antenna array layout definition from a YAML (.yaml, .yml) or JSON (.json) file.

    Parameters
    ----------
    layout_path : str
        Path to the layout YAML or JSON file.

    Returns
    -------
    dict or list of dict
        Layout dictionary containing 'metadata' and 'elements', or list of element dictionaries.
    """
    layout_path = os.path.abspath(layout_path)
    ext = os.path.splitext(layout_path)[1].lower()

    with open(layout_path, "r", encoding="utf-8") as f:
        if ext in [".yaml", ".yml"]:
            import yaml

            return yaml.safe_load(f)
        elif ext == ".json":
            import json

            return json.load(f)
        else:
            try:
                import yaml

                return yaml.safe_load(f)
            except Exception:
                f.seek(0)
                import json

                return json.load(f)


def main(argv=None):
    """Command-line entry point (installed as `hfss_array_builder`)."""
    import argparse

    parser = argparse.ArgumentParser(description="HFSS Full-Wave Array Synthesis CLI.")
    parser.add_argument("project_path", help="Path to the AEDT project file (.aedt)")
    parser.add_argument("source_design_name", help="Name of the unit-cell source design")
    parser.add_argument(
        "layout_path",
        nargs="?",
        default=None,
        help="Optional path to custom elements layout YAML or JSON file",
    )
    parser.add_argument(
        "-f",
        "--overwrite",
        "--overwrite-design",
        action="store_true",
        default=False,
        help="Automatically overwrite/clear target HFSS design if it already exists",
    )
    parser.add_argument(
        "-y",
        "--yes",
        "--accept-all",
        action="store_true",
        default=False,
        help="Accept all affirmative defaults (use existing CS, overwrite existing design, setup results, and run simulation)",
    )
    parser.add_argument(
        "--build-only",
        action="store_true",
        default=False,
        help="Build model only: automatically decline results creation and decline simulation run without prompting",
    )
    parser.add_argument(
        "--results-only",
        action="store_true",
        default=False,
        help="Skip building; only set up post-processing results on an existing design",
    )
    parser.add_argument(
        "--simulate-only",
        action="store_true",
        default=False,
        help="Skip building; only launch simulation ('Analyze All') on an existing design",
    )
    parser.add_argument(
        "--use-existing-cs",
        "--use-existing-phase-centre",
        action="store_true",
        default=False,
        help="Use existing PhaseCentreCS in the unit cell without prompting or re-fitting it",
    )
    parser.add_argument(
        "--centre-freq",
        "--center-freq",
        type=float,
        default=None,
        help=(
            "Centre frequency in GHz (default: layout metadata, else the unit-cell/array setup, "
            f"else {DEFAULT_CENTRE_FREQUENCY_GHZ})"
        ),
    )
    parser.add_argument(
        "--bandwidth", type=float, default=4.0, help="Sweep bandwidth in GHz (default: 4.0)"
    )
    parser.add_argument(
        "--phase-reference-mode",
        "--phase-ref-mode",
        choices=["single", "per-port"],
        default="single",
        help="Far-field phase reference mode for SBR+ antenna linking: 'single' (array centroid PhaseCentreCS, default) or 'per-port' (individual element CS)",
    )

    parser.add_argument(
        "--port",
        type=int,
        default=None,
        help=(
            "gRPC port of the AEDT session to connect to; a new session is started on it if "
            "none answers (default: the running AEDT GUI's port, else 50051)"
        ),
    )

    args = parser.parse_args(argv)

    # Determine target design name from layout file if provided
    target_design_name = None
    if args.layout_path:
        base_name = os.path.splitext(os.path.basename(args.layout_path))[0]
        words = base_name.replace("-", "_").split("_")
        layout_suffix = "".join(w.capitalize() for w in words if w)
        target_design_name = f"{args.source_design_name}{layout_suffix}"

    builder = MimoHfssBuilder(
        project_path=args.project_path,
        source_design_name=args.source_design_name,
        target_design_name=target_design_name,
        centre_frequency_ghz=args.centre_freq,
        bandwidth_ghz=args.bandwidth,
        phase_reference_mode=args.phase_reference_mode,
        grpc_port=args.port,
        non_graphical=True,
    )

    try:
        builder.connect_desktop()

        if args.results_only or args.simulate_only:
            target_design = builder.target_design_name
            print(
                f"\n[INFO] Connecting directly to target design '{target_design}' in '{args.project_path}'..."
            )
            builder.target_design_app = Hfss(
                project=args.project_path,
                design=target_design,
                new_desktop=False,
                close_on_exit=False,
            )

            if args.results_only:
                builder.create_post_processing_reports()
                builder.target_design_app.save_project()

            if args.simulate_only:
                print("\nRunning HFSS built-in design validation...")
                validation_ok = builder.target_design_app.validate_simple()
                if validation_ok == 1 or validation_ok is True:
                    print("  Design validation: PASSED.")
                else:
                    print("  [WARNING] Design validation returned issues.")
                print("\nLaunching HFSS simulation ('Analyze All')...")
                builder.target_design_app.analyze()
                builder.target_design_app.save_project()
        else:
            if args.layout_path:
                print(f"Loading custom layout from: {args.layout_path}")
                elements_list = load_layout_file(args.layout_path)
            else:
                layout_freq = args.centre_freq or DEFAULT_CENTRE_FREQUENCY_GHZ
                print(
                    f"No custom layout provided. Using default coplanar layout at {layout_freq} GHz..."
                )
                elements_list = builder.calculate_default_coplanar_layout(
                    transmitter_count=4,
                    receiver_count=4,
                    operating_frequency_ghz=layout_freq,
                    subarray_spacing_mm=10.0,
                )

            if args.build_only:
                setup_results = False
                run_simulation = False
            elif args.yes:
                setup_results = True
                run_simulation = True
            else:
                setup_results = None
                run_simulation = None

            builder.synthesize_array_in_hfss(
                elements_list,
                use_existing_cs=True if (args.use_existing_cs or args.yes) else None,
                overwrite=True if (args.overwrite or args.yes) else None,
                setup_results=setup_results,
                run_simulation=run_simulation,
            )
    except BuildAbortedError as e:
        print(f"[INFO] {e}")
        return 0
    except NativeAedtRefusedError as e:
        print(f"[ERROR] {e}")
        return 2
    finally:
        builder.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
