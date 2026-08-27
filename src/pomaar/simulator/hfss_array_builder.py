#!/usr/bin/env python3
"""
Module 1: HFSS Full-Wave Array Synthesis.
Programmatically builds a planar MIMO antenna array layout on a single PCB substrate in HFSS
using a template-based, boolean-driven assembly workflow.
"""

import math
import os
import sys

from ansys.aedt.core import Desktop, Hfss
from ansys.aedt.core.generic.constants import Axis, Gravity


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
        grpc_port=50051,
        non_graphical=True,
        centre_frequency_ghz=79.0,
        bandwidth_ghz=4.0,
    ):
        self.project_path = os.path.abspath(project_path)
        self.source_design_name = source_design_name
        self.target_design_name = target_design_name if target_design_name else f"{source_design_name}MimoArray"
        self.pcb_margin_mm = pcb_margin_mm
        self.grpc_port = grpc_port
        self.non_graphical = non_graphical
        self.centre_frequency_ghz = centre_frequency_ghz
        self.bandwidth_ghz = bandwidth_ghz

        self.desktop_session = None
        self.source_design_app = None
        self.target_design_app = None
        self.is_new_desktop = False
        self.is_new_design = True

    def connect_desktop(self):
        """Starts or connects to an existing AEDT session via gRPC."""
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
        metric_choice=None,
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
        """
        if isinstance(elements, dict):
            layout_metadata = elements.get("metadata", {})
            elements_list = elements.get("elements", [])
        else:
            layout_metadata = {}
            elements_list = elements

        if operating_frequency_ghz is not None:
            self.centre_frequency_ghz = operating_frequency_ghz
        elif "center_frequency_ghz" in layout_metadata:
            self.centre_frequency_ghz = float(layout_metadata["center_frequency_ghz"])
            operating_frequency_ghz = self.centre_frequency_ghz
        else:
            operating_frequency_ghz = self.centre_frequency_ghz

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
            if obj.material_name.lower() == "vacuum" or "RadiatingSurface" in name or "Airbox" in name or "RadiationBox" in name:
                print(f"  [Skipped Boundary] {name} ({obj.material_name})")
                continue

            if name.startswith("L12_") or name.startswith("L23_") or name.startswith("L34_") or name.startswith("L45_"):
                bbox = obj.bounding_box
                z_coords = [float(bbox[2]), float(bbox[5])]
                global_layers[name] = {
                    "material": obj.material_name,
                    "z_min": min(z_coords),
                    "z_max": max(z_coords),
                }
                print(f"  [Global Layer]     {name} ({obj.material_name}, Z=[{min(z_coords):.2f}, {max(z_coords):.2f}] mm)")

            elif name.startswith("L2_Ground") or name.startswith("Ground_Plane") or name.startswith("GND"):
                bbox = obj.bounding_box
                z_coords = [float(bbox[2]), float(bbox[5])]
                global_layers[name] = {
                    "material": obj.material_name,
                    "z_min": min(z_coords),
                    "z_max": max(z_coords),
                }
                print(f"  [Global Ground]    {name} ({obj.material_name}, Z=[{min(z_coords):.2f}, {max(z_coords):.2f}] mm)")

            elif "PortSheet" in name or name.startswith("Port_"):
                port_sheets.append(name)
                print(f"  [Port Sheet]       {name}")

            else:
                active_elements.append(name)
                print(f"  [Active Element]   {name}")

        run_opt = False
        offset = [0.0, 0.0, 0.0]

        cs_map = {cs.name: cs for cs in source_modeler.coordinate_systems}
        if "PhaseCentreCS" in cs_map:
            cs_obj = cs_map["PhaseCentreCS"]
            try:
                dx_mm = float(self.source_design_app.evaluate_expression(cs_obj.props["OriginX"])) * 1000.0
                dy_mm = float(self.source_design_app.evaluate_expression(cs_obj.props["OriginY"])) * 1000.0
                dz_mm = float(self.source_design_app.evaluate_expression(cs_obj.props["OriginZ"])) * 1000.0
                offset = [round(dx_mm, 3), round(dy_mm, 3), round(dz_mm, 3)]
            except Exception:
                try:
                    offset = [float(val) for val in cs_obj.origin]
                except Exception:
                    offset = [0.0, 0.0, 0.0]
            print(f"\n[INFO] PhaseCentreCS found in unit-cell design. Origin offset: {offset} mm.")
            if use_existing_cs is True:
                user_input = "y"
            elif use_existing_cs is False:
                user_input = "n"
            elif sys.stdin.isatty():
                user_input = input("Do you want to use the existing PhaseCentreCS? (y/n) [default: y]: ")
            else:
                user_input = "y"
            if user_input.strip().lower() in ["n", "no"]:
                print("[INFO] Re-calculating Phase Centre using Optimetrics...")
                run_opt = True
        else:
            print("\n" + "=" * 80)
            print("[WARNING] PhaseCentreCS was NOT found in the unit-cell design!")
            print("=" * 80 + "\n")
            if use_existing_cs is True:
                run_opt = False
                print("[INFO] Proceeding without PhaseCentreCS offset [0.0, 0.0, 0.0] mm.")
            elif sys.stdin.isatty():
                user_input = input("Do you want to run Optimetrics to calculate the Phase Centre? (y/n) [default: y]: ")
                if user_input.strip().lower() not in ["n", "no"]:
                    run_opt = True
                else:
                    user_input2 = input("Proceed without PhaseCentreCS offset? (y/n) [default: y]: ")
                    if user_input2.strip().lower() in ["n", "no"]:
                        print("[INFO] Aborted by user.")
                        self.close()
                        sys.exit(0)
                    print("[INFO] Proceeding with zero offset [0.0, 0.0, 0.0] mm.")
            else:
                run_opt = True

        if run_opt:
            offset = self._run_phase_centre_opt(source_modeler, global_layers, active_elements, operating_frequency_ghz)

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
            if overwrite is True:
                user_input = "y"
            elif overwrite is False:
                user_input = "n"
            elif sys.stdin.isatty():
                user_input = input(
                    (
                        f"Design '{self.target_design_name}' already exists. "
                        "Do you want to clear it and overwrite? (y/n) [default: y]: "
                    )
                )
            else:
                user_input = "y"
            if user_input.strip().lower() in ["n", "no"]:
                print("[INFO] Aborted by user to prevent overwriting existing design.")
                self.close()
                sys.exit(0)

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
            for var_name, var_expr in self.source_design_app.variable_manager.design_variables.items():
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
                    uc_x_extents.extend([abs(float(bb[0]) - offset[0]), abs(float(bb[3]) - offset[0])])
                    uc_y_extents.extend([abs(float(bb[1]) - offset[1]), abs(float(bb[4]) - offset[1])])
            except Exception:
                pass

        unit_cell_extent_x = max(uc_x_extents) if uc_x_extents else 5.0
        unit_cell_extent_y = max(uc_y_extents) if uc_y_extents else 5.0
        unit_cell_extent_max = max(unit_cell_extent_x, unit_cell_extent_y)

        # Calculate clearance using lambda_0 / 4 (quarter wavelength) rule of thumb
        speed_of_light_mm_s = 2.99792458e11
        wavelength_mm = speed_of_light_mm_s / (operating_frequency_ghz * 1e9)
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

        for layer_name in sorted_layers:
            layer_info = global_layers[layer_name]
            layer_obj = self.source_design_app.modeler.get_object_from_name(layer_name)
            bbox = layer_obj.bounding_box
            x_min_uc, x_max_uc = float(bbox[0]), float(bbox[3])

            y_min_active_uc = None
            y_max_active_uc = None
            for active_name in active_elements + port_sheets:
                try:
                    act_obj = self.source_design_app.modeler.get_object_from_name(active_name)
                    act_bbox = act_obj.bounding_box
                    act_y_min, act_y_max = float(act_bbox[1]), float(act_bbox[4])
                    if y_min_active_uc is None or act_y_min < y_min_active_uc:
                        y_min_active_uc = act_y_min
                    if y_max_active_uc is None or act_y_max > y_max_active_uc:
                        y_max_active_uc = act_y_max
                except Exception:
                    pass

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
                yaw_deg = float(element.get("yaw", element.get("rotation_yaw", element.get("rotation", 0.0))))
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
                print(f"  Creating substrate '{layer_name}' (material={material}, thickness={thickness:.2f} mm)")
                sub_board = target_modeler.create_box(
                    origin=["-arrayBoardWidth / 2", "-arrayBoardLength / 2", f"{z_min - offset[2]:.4f}mm"],
                    sizes=["arrayBoardWidth", "arrayBoardLength", f"{thickness:.4f}mm"],
                    name=layer_name,
                    material=material,
                )
                sub_board.transparency = 0.5

            elif layer_name.endswith("_Ground"):
                ground_layer_names.append(layer_name)
                if thickness == 0.0:
                    print(f"  Creating ground plane sheet '{layer_name}' at Z={z_min - offset[2]:.2f} mm")
                    target_modeler.create_rectangle(
                        orientation="XY",
                        origin=["-arrayBoardWidth / 2", "-arrayBoardLength / 2", f"{z_min - offset[2]:.4f}mm"],
                        sizes=["arrayBoardWidth", "arrayBoardLength"],
                        name=layer_name,
                    )
                else:
                    print(f"  Creating ground plane block '{layer_name}' (thickness={thickness:.2f} mm)")
                    target_modeler.create_box(
                        origin=["-arrayBoardWidth / 2", "-arrayBoardLength / 2", f"{z_min - offset[2]:.4f}mm"],
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
            element_yaw = float(element.get("yaw", element.get("rotation_yaw", element.get("rotation", 0.0))))
            pol = str(element.get("polarization", element.get("pol", ""))).strip().upper()

            # Append polarization suffix if specified and not already present in the label
            if pol and not raw_label.upper().endswith(f"_{pol}"):
                label = f"{raw_label}_{pol}"
            else:
                label = raw_label

            cs_name = f"CS_{label}"
            print(f"\n  Setting up element {label} (yaw={element_yaw:.1f} deg)...")

            # Determine parametric origin of the element CS
            if pos_expr:
                cs_origin = [str(pos_expr[0]), str(pos_expr[1]), f"{-offset[2]:.4f}mm"]
            else:
                cs_origin = [f"{pos[0]:.4f}mm", f"{pos[1]:.4f}mm", f"{-offset[2]:.4f}mm"]

            # Calculate pointing vectors for the CS rotation around Z
            rad = math.radians(element_yaw)
            cos_val = math.cos(rad)
            sin_val = math.sin(rad)
            x_pointing = [cos_val, sin_val, 0.0]
            y_pointing = [-sin_val, cos_val, 0.0]

            # Recreate coordinate system if it exists
            cs_map = {cs.name: cs for cs in target_modeler.coordinate_systems}
            if cs_name in cs_map:
                try:
                    cs_map[cs_name].delete()
                except Exception:
                    pass

            print(f"  Creating relative coordinate system '{cs_name}' at origin {cs_origin}...")
            target_modeler.create_coordinate_system(
                origin=cs_origin,
                reference_cs="Global",
                name=cs_name,
                mode="axis",
                x_pointing=x_pointing,
                y_pointing=y_pointing,
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

                # Track dummy solids
                for key, dummy_src_list in dummy_solids.items():
                    if base_name in dummy_src_list:
                        replicated_dummy_mapping[key].append(new_name)

            # Rotate element if needed (centred at the origin)
            if element_yaw != 0.0:
                target_modeler.rotate(
                    assignment=renamed_objs,
                    axis=Axis.Z,
                    angle=element_yaw,
                )

            # Calculate the rotated phase centre offset
            rot_dx = offset[0] * cos_val - offset[1] * sin_val
            rot_dy = offset[0] * sin_val + offset[1] * cos_val
            rot_dz = offset[2]

            # In Z, the element must be translated by (pos[2] - rot_dz) so its vertical positioning
            # relative to the substrate and ground plane matches the unit-cell design exactly.
            z_pos_val = float(pos[2]) if len(pos) > 2 else 0.0
            z_shift_str = f"{z_pos_val - rot_dz:.4f}mm"

            # Move element to final position
            if pos_expr:
                x_move = f"{pos_expr[0]} - {rot_dx:.4f}mm" if abs(rot_dx) > 1e-4 else str(pos_expr[0])
                y_move = f"{pos_expr[1]} - {rot_dy:.4f}mm" if abs(rot_dy) > 1e-4 else str(pos_expr[1])
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
                print(f"  Warning: Target solid '{target_solid}' not found in layout. Skipping boolean {operation}.")
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

        for port_sheet in replicated_ports_list:
            suffix = port_sheet.replace("PortSheet", "")
            if suffix.startswith("_"):
                suffix = suffix[1:]
            port_name = f"Port_{suffix}"

            print(f"  Assigning wave port to sheet '{port_sheet}' referencing '{default_ground}'")
            try:
                self.target_design_app.wave_port(
                    assignment=port_sheet,
                    reference=default_ground,
                    integration_line=Gravity.ZNeg,
                    name=port_name,
                )
            except Exception as e:
                print(f"  Warning: wave_port assignment failed ({e}), falling back to lumped_port...")
                self.target_design_app.lumped_port(
                    assignment=port_sheet,
                    reference=default_ground,
                    integration_line=Gravity.ZNeg,
                    impedance=50.0,
                    name=port_name,
                )

        # Enable Port Post Processing Effects to satisfy warning constraints
        try:
            port_assignments = {}
            for port_sheet in replicated_ports_list:
                suffix = port_sheet.replace("PortSheet", "")
                if suffix.startswith("_"):
                    suffix = suffix[1:]
                port_name = f"Port_{suffix}"
                port_assignments[f"{port_name}:1"] = ("1W", "0deg")

            self.target_design_app.edit_sources(assignment=port_assignments, include_port_post_processing=True)
            print("  Enabled 'Include Port Post Processing Effects' in Edit Sources.")
        except Exception as e:
            print(f"  Warning: Could not enable port post processing effects ({e})")

        # --- STEP 7: Create Radiation Airbox ---
        print("\nCreating radiation boundary Airbox...")
        # Get overall board bounds in Z direction
        all_z_coords = []
        for layer_info in global_layers.values():
            all_z_coords.extend([layer_info["z_min"], layer_info["z_max"]])
        overall_z_min = (min(all_z_coords) if all_z_coords else 0.0) - offset[2]
        overall_z_max = (max(all_z_coords) if all_z_coords else 1.0) - offset[2]
        total_z_thickness = overall_z_max - overall_z_min

        print(f"  Creating flush lateral Airbox solid (airboxClearance={airbox_clearance_mm:.2f} mm)")
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
        airbox_obj.transparency = 0.9
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

        # Configure simulation setup and sweep
        if self.is_new_design:
            # Configure simulation setup and sweep in the target design
            self.configure_simulation_setup()
        else:
            print(
                "\n[INFO] Reusing existing design: Preserving all existing simulation setups and frequency sweeps."
            )

        # As the last step of building, ask the user whether to set up results (far-field sphere & post-processing reports)
        if setup_results is None:
            if sys.stdin.isatty():
                user_input = input(
                    "\nDo you want to set up results (far-field sphere & post-processing reports)? (y/n) [default: y]: "
                )
                setup_results = user_input.strip().lower() not in ["n", "no"]
            else:
                setup_results = True

        if setup_results:
            self.create_post_processing_reports(metric_choice=metric_choice)
        else:
            print("\n[INFO] Skipping results setup as requested.")

        # Run Design Validation
        print("\nRunning HFSS built-in design validation...")
        validation_ok = self.target_design_app.validate_simple()
        if validation_ok == 1 or validation_ok is True:
            print("  Design validation: PASSED.")
        else:
            print("  [ERROR] Design validation: FAILED!")
            assert False, "HFSS design validation failed."

        self.target_design_app.save_project()
        print(f"\nMIMO array synthesis successfully completed in design '{self.target_design_name}'.")

        # As the very last step, ask the user if they want to launch the simulation ('Analyze All')
        if run_simulation is None:
            if sys.stdin.isatty():
                user_input = input(
                    "\nDo you want to launch the simulation ('Analyze All')? (y/n) [default: y]: "
                )
                run_simulation = user_input.strip().lower() not in ["n", "no"]
            else:
                run_simulation = True

        if run_simulation:
            print("\nLaunching HFSS simulation ('Analyze All')...")
            self.target_design_app.analyze()
        else:
            print("\n[INFO] Skipping simulation execution. Model synthesis is complete.")

    def configure_simulation_setup(self):
        """Configures the single-frequency adaptive mesh setup and frequency sweep in the target design."""
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

        # Check and purge any existing analysis setups and their sweeps
        # Check and purge any existing analysis setups and their sweeps
        try:
            setup_names = list(self.target_design_app.setup_names)
        except Exception:
            setup_names = []

        if setup_names:
            print(f"  Clearing {len(setup_names)} existing analysis setups and their sweeps...")
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

        setup_name = "ArraySetup"
        print(
            f"\nConfiguring single-frequency adaptive mesh setup '{setup_name}' at {self.centre_frequency_ghz} GHz..."
        )

        # Create setup fresh
        setup = self.target_design_app.create_setup(name=setup_name)

        setup.props["Frequency"] = f"{self.centre_frequency_ghz}GHz"
        setup.props["MaximumPasses"] = 21
        setup.props["MaxDeltaS"] = 0.02
        setup.props["SaveFields"] = True
        setup.update()

        # Deduce frequency sweep bounds
        start_freq = self.centre_frequency_ghz - self.bandwidth_ghz / 2.0
        end_freq = self.centre_frequency_ghz + self.bandwidth_ghz / 2.0

        sweep_name = "Sweep"
        print(
            f"Configuring interpolating frequency sweep '{sweep_name}' ({start_freq:.2f} GHz - {end_freq:.2f} GHz, 401",
            "points)...",
        )

        # Try to delete default sweep if created automatically
        try:
            # Initialize sweeps cache to prevent NoneType iterable error
            _ = setup.sweeps
            for sname in list(setup.get_sweep_names()):
                setup.delete_sweep(sname)
        except Exception:
            pass

        self.target_design_app.create_linear_count_sweep(
            setup=setup_name,
            unit="GHz",
            start_frequency=start_freq,
            stop_frequency=end_freq,
            num_of_freq_points=401,
            name=sweep_name,
            save_fields=True,
            sweep_type="Interpolating",
        )

    def create_post_processing_reports(self, metric_choice=None):
        """Creates standard S-parameter and Far-Field reports in the target design."""
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

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
                print(f"  Generating Tx1-to-Rx crosstalk report...")
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

        # --- STEP B: Prompt user for Radiation Pattern Metric ---
        if metric_choice is None:
            if sys.stdin.isatty():
                print("\nWhich radiation pattern metric to use?")
                print("  1. Directivity")
                print("  2. Gain")
                print("  3. Realized gain")
                user_metric = input("Select option (1/2/3) [default: 3]: ")
                metric_choice = user_metric.strip()
            else:
                metric_choice = "3"

        if metric_choice == "1" or str(metric_choice).lower() in ["directivity", "1"]:
            metric_name = "Directivity"
            total_qty = "DirTotal"
            copolar_qty = "DirCoPolar"
            crosspolar_qty = "DirCrossPolar"
        elif metric_choice == "2" or str(metric_choice).lower() in ["gain", "2"]:
            metric_name = "Gain"
            total_qty = "GainTotal"
            copolar_qty = "GainCoPolar"
            crosspolar_qty = "GainCrossPolar"
        else:
            metric_name = "Realized gain"
            total_qty = "RealizedGainTotal"
            copolar_qty = "RealizedGainCoPolar"
            crosspolar_qty = "RealizedGainCrossPolar"

        print(f"  Using radiation pattern metric: '{metric_name}' ({total_qty})")

        # 5. Far Field Setups and Radiation Patterns
        # (a) InfiniteSphere: Standard IEEE Theta-Phi system (z-axis zenith, theta 0..180, phi 0..360)
        sphere_name = "InfiniteSphere"
        print(f"  Configuring Infinite Sphere '{sphere_name}' (IEEE convention: theta 0..180 deg, phi 0..360 deg)...")
        try:
            self.target_design_app.insert_infinite_sphere(
                name=sphere_name,
                definition="Theta-Phi",
                phi_start=0, phi_stop=360, phi_step=5,
                theta_start=0, theta_stop=180, theta_step=1,
                units="deg"
            )
        except Exception as e:
            print(f"  Warning: Failed to create infinite sphere '{sphere_name}' ({e})")

        # (b) Setup 1: 'Phi0 cut' in Az Over El system (azimuth -180..180 step 2, elevation 0..0 step 0)
        phi0_sphere = "Phi0 cut"
        print(f"  Configuring Infinite Sphere '{phi0_sphere}' (Az Over El: az -180..180 step 2, el 0..0 step 0)...")
        try:
            self.target_design_app.insert_infinite_sphere(
                name=phi0_sphere,
                definition="Az Over El",
                phi_start=-180, phi_stop=180, phi_step=2,
                theta_start=0, theta_stop=0, theta_step=0,
                units="deg"
            )
        except Exception as e:
            print(f"  Warning: Failed to create infinite sphere '{phi0_sphere}' ({e})")

        # (b) Setup 2: 'Phi90 cut' in Az Over El system (azimuth 0..0 step 0, elevation -180..180 step 2)
        phi90_sphere = "Phi90 cut"
        print(f"  Configuring Infinite Sphere '{phi90_sphere}' (Az Over El: az 0..0 step 0, el -180..180 step 2)...")
        try:
            self.target_design_app.insert_infinite_sphere(
                name=phi90_sphere,
                definition="Az Over El",
                phi_start=0, phi_stop=0, phi_step=0,
                theta_start=-180, theta_stop=180, theta_step=2,
                units="deg"
            )
        except Exception as e:
            print(f"  Warning: Failed to create infinite sphere '{phi90_sphere}' ({e})")

        # Create radiation pattern reports
        # --- 3D Pattern (3D Polar Plot on InfiniteSphere) ---
        plot_name = f"{metric_name} 3D"
        if plot_name not in existing_reports:
            try:
                print(f"  Generating {plot_name} report...")
                vars_3d = {"Theta": ["All"], "Phi": ["All"], "Freq": [f"{self.centre_frequency_ghz}GHz"]}
                self.target_design_app.post.create_report(
                    expressions=[f"db({total_qty})"],
                    setup_sweep_name=setup_sweep,
                    variations=vars_3d,
                    primary_sweep_variable="Phi",
                    secondary_sweep_variable="Theta",
                    report_category="Far Fields",
                    plot_name=plot_name,
                    context=sphere_name,
                    plot_type="3D Polar Plot"
                )
            except Exception as e:
                print(f"  Warning: Failed to create {plot_name} report ({e})")
        else:
            print(f"  Preserving existing report '{plot_name}'")

        # (c) Rectangular Far-Field Cut 1: Phi0 geometry with primary sweep 'AzimuthAngle'
        plot_name = f"{metric_name} Phi0"
        if plot_name not in existing_reports:
            try:
                print(f"  Generating {plot_name} report (Az Over El)...")
                vars_phi0 = {"AzimuthAngle": ["All"], "ElevationAngle": ["All"], "Freq": [f"{self.centre_frequency_ghz}GHz"]}
                self.target_design_app.post.create_report(
                    expressions=[f"db({total_qty})"],
                    setup_sweep_name=setup_sweep,
                    variations=vars_phi0,
                    primary_sweep_variable="AzimuthAngle",
                    report_category="Far Fields",
                    plot_name=plot_name,
                    context=phi0_sphere,
                    plot_type="Rectangular Plot"
                )
            except Exception as e:
                print(f"  Warning: Failed to create {plot_name} report ({e})")
        else:
            print(f"  Preserving existing report '{plot_name}'")

        # (c) Rectangular Far-Field Cut 2: Phi90 geometry with primary sweep 'ElevationAngle'
        plot_name = f"{metric_name} Phi90"
        if plot_name not in existing_reports:
            try:
                print(f"  Generating {plot_name} report (Az Over El)...")
                vars_phi90 = {"AzimuthAngle": ["All"], "ElevationAngle": ["All"], "Freq": [f"{self.centre_frequency_ghz}GHz"]}
                self.target_design_app.post.create_report(
                    expressions=[f"db({total_qty})"],
                    setup_sweep_name=setup_sweep,
                    variations=vars_phi90,
                    primary_sweep_variable="ElevationAngle",
                    report_category="Far Fields",
                    plot_name=plot_name,
                    context=phi90_sphere,
                    plot_type="Rectangular Plot"
                )
            except Exception as e:
                print(f"  Warning: Failed to create {plot_name} report ({e})")
        else:
            print(f"  Preserving existing report '{plot_name}'")

        # (d) XPD metric defined as CoPolar - CrossPolar (Ludwig-3 definition)
        xpd_expression = f"db({copolar_qty}) - db({crosspolar_qty})"

        plot_name = "XPD Phi0"
        if plot_name not in existing_reports:
            try:
                print(f"  Generating XPD Phi0 report...")
                vars_phi0 = {"AzimuthAngle": ["All"], "ElevationAngle": ["All"], "Freq": [f"{self.centre_frequency_ghz}GHz"]}
                self.target_design_app.post.create_report(
                    expressions=xpd_expression,
                    setup_sweep_name=setup_sweep,
                    variations=vars_phi0,
                    primary_sweep_variable="AzimuthAngle",
                    report_category="Far Fields",
                    plot_name=plot_name,
                    context=phi0_sphere,
                    plot_type="Rectangular Plot"
                )
            except Exception as e:
                print(f"  Warning: Failed to create {plot_name} report ({e})")
        else:
            print(f"  Preserving existing report '{plot_name}'")

        plot_name = "XPD Phi90"
        if plot_name not in existing_reports:
            try:
                print(f"  Generating XPD Phi90 report...")
                vars_phi90 = {"AzimuthAngle": ["All"], "ElevationAngle": ["All"], "Freq": [f"{self.centre_frequency_ghz}GHz"]}
                self.target_design_app.post.create_report(
                    expressions=xpd_expression,
                    setup_sweep_name=setup_sweep,
                    variations=vars_phi90,
                    primary_sweep_variable="ElevationAngle",
                    report_category="Far Fields",
                    plot_name=plot_name,
                    context=phi90_sphere,
                    plot_type="Rectangular Plot"
                )
            except Exception as e:
                print(f"  Warning: Failed to create {plot_name} report ({e})")
        else:
            print(f"  Preserving existing report '{plot_name}'")

    def run_solve(self):
        """Triggers the HFSS simulation setup and solve using 'Analyze All'."""
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

        # Ensure the setup is configured
        self.configure_simulation_setup()

        print("Solving full-wave HFSS array design using 'Analyze All'...")
        self.target_design_app.analyze()

    def export_coupling_s_parameters(self, output_touchstone_path):
        """Exports solved S-parameter coupling matrix to a Touchstone (.sNp) file."""
        if not self.target_design_app:
            raise RuntimeError("Array design is not resolved.")

        output_touchstone_path = os.path.abspath(output_touchstone_path)
        print(f"Exporting coupling S-parameters to: {output_touchstone_path}")

        self.target_design_app.export_touchstone(
            setup_name="ArraySetup", sweep_name="LastSweep", filename=output_touchstone_path
        )

    def _run_phase_centre_opt(self, source_modeler, global_layers, active_elements, operating_frequency_ghz):
        """Runs HFSS Optimetrics to find the exact phase centre of the antenna."""
        print("\nStarting automated Phase Centre extraction...")

        # Calculate substrate span and centre
        x_min_uc, x_max_uc = None, None
        y_min_uc, y_max_uc = None, None
        for obj_name, layer_info in global_layers.items():
            if "_substrate" in obj_name.lower():
                layer_obj = source_modeler[obj_name]
                bbox = layer_obj.bounding_box
                cx_min = float(bbox[0])
                cx_max = float(bbox[3])
                cy_min = float(bbox[1])
                cy_max = float(bbox[4])
                if x_min_uc is None or cx_min < x_min_uc:
                    x_min_uc = cx_min
                if x_max_uc is None or cx_max > x_max_uc:
                    x_max_uc = cx_max
                if y_min_uc is None or cy_min < y_min_uc:
                    y_min_uc = cy_min
                if y_max_uc is None or cy_max > y_max_uc:
                    y_max_uc = cy_max

        # Fallback values
        if x_min_uc is None:
            x_min_uc, x_max_uc = -5.0, 5.0
        if y_min_uc is None:
            y_min_uc, y_max_uc = -5.0, 5.0

        x_span = x_max_uc - x_min_uc
        y_span = y_max_uc - y_min_uc
        x_centre = (x_min_uc + x_max_uc) / 2.0
        y_centre = (y_min_uc + y_max_uc) / 2.0

        # Find top Z coordinate
        top_z = 0.0
        for layer_info in global_layers.values():
            if layer_info["z_max"] > top_z:
                top_z = layer_info["z_max"]
        for obj_name in active_elements:
            layer_obj = source_modeler[obj_name]
            bbox = layer_obj.bounding_box
            z_max = float(bbox[5])
            if z_max > top_z:
                top_z = z_max

        print(f"  Unit-cell lateral spans: X={x_span:.2f} mm, Y={y_span:.2f} mm")
        print(f"  Top-most Z layer: {top_z:.3f} mm")

        # Set up design variables for optimization
        print("  Creating design variables...")
        self.source_design_app["PhaseCentreX"] = f"{x_centre:.3f}mm"
        self.source_design_app["PhaseCentreY"] = f"{y_centre:.3f}mm"
        self.source_design_app["PhaseCentreZ"] = f"{top_z:.3f}mm"

        # Activate variable optimization with ranges
        self.source_design_app.activate_variable_optimization(
            "PhaseCentreX", minimum=f"{x_centre - x_span / 2.0:.3f}mm", maximum=f"{x_centre + x_span / 2.0:.3f}mm"
        )
        self.source_design_app.activate_variable_optimization(
            "PhaseCentreY", minimum=f"{y_centre - y_span / 2.0:.3f}mm", maximum=f"{y_centre + y_span / 2.0:.3f}mm"
        )
        self.source_design_app.activate_variable_optimization(
            "PhaseCentreZ", minimum=f"{top_z - 1.0:.3f}mm", maximum=f"{top_z + 2.0:.3f}mm"
        )

        # Create Relative Coordinate System
        temp_cs_name = "PhaseCentreCS_Opt"
        cs_map = {cs.name: cs for cs in source_modeler.coordinate_systems}
        if temp_cs_name in cs_map:
            try:
                cs_map[temp_cs_name].delete()
            except Exception:
                pass

        print(f"  Creating temporary coordinate system '{temp_cs_name}'...")
        source_modeler.create_coordinate_system(
            origin=["PhaseCentreX", "PhaseCentreY", "PhaseCentreZ"], reference_cs="Global", name=temp_cs_name
        )

        # Create Far Field Infinite Sphere Setup
        # We sweep theta from -40 to 40 degrees at phi=0
        sphere_name = "PhaseCentreSphere"
        print(f"  Creating far-field infinite sphere setup '{sphere_name}' bound to '{temp_cs_name}'...")
        try:
            self.source_design_app.field_setups[sphere_name].delete()
        except Exception:
            pass

        self.source_design_app.insert_infinite_sphere(
            name=sphere_name,
            phi_start=0,
            phi_stop=0,
            phi_step=1,
            theta_start=-40,
            theta_stop=40,
            theta_step=2,
            units="deg",
            custom_coordinate_system=temp_cs_name,
        )

        # Create Optimetrics Optimization Setup
        opt_setup_name = "PhaseCentreOpt"
        print(f"  Creating Optimetrics setup '{opt_setup_name}'...")
        try:
            self.source_design_app.optimizations.delete(opt_setup_name)
        except Exception:
            pass

        source_setup = self.source_design_app.setup_names[0] if self.source_design_app.setup_names else None
        _ = self.source_design_app.optimizations.add(
            calculation="pk2pk(cang_deg(rEphi))",
            ranges={
                "Theta": ("-40deg", "40deg"),
                "Phi": "0deg",
                "Freq": f"{operating_frequency_ghz}GHz",
            },
            optimization_type="Optimization",
            variables=["PhaseCentreX", "PhaseCentreY", "PhaseCentreZ"],
            name=opt_setup_name,
            context=sphere_name,
            report_type="Far Fields",
            condition="Minimize",
            solution=source_setup,
        )

        print("  Running Phase Centre Optimization in HFSS...")
        try:
            self.source_design_app.analyze_setup(opt_setup_name)
        except Exception as e:
            print(f"  [ERROR] Optimization failed: {e}")
            print("  Falling back to zero offset [0.0, 0.0, 0.0] mm.")
            return [0.0, 0.0, 0.0]

        # Retrieve optimized values
        print("  Retrieving optimized coordinates...")
        # evaluate_expression returns values in SI units (meters). Convert to millimeters.
        opt_x = float(self.source_design_app.evaluate_expression("PhaseCentreX")) * 1000.0
        opt_y = float(self.source_design_app.evaluate_expression("PhaseCentreY")) * 1000.0
        opt_z = float(self.source_design_app.evaluate_expression("PhaseCentreZ")) * 1000.0

        # Round to 3 decimal places
        opt_x_rounded = round(opt_x, 3)
        opt_y_rounded = round(opt_y, 3)
        opt_z_rounded = round(opt_z, 3)

        print(f"  Optimized Phase Centre coordinates: [{opt_x_rounded}, {opt_y_rounded}, {opt_z_rounded}] mm")

        # Clean up temporary setups
        print("  Cleaning up temporary Optimetrics configurations...")
        try:
            self.source_design_app.optimizations.delete(opt_setup_name)
        except Exception:
            pass
        try:
            self.source_design_app.field_setups[sphere_name].delete()
        except Exception:
            pass
        try:
            cs_map = {cs.name: cs for cs in source_modeler.coordinate_systems}
            if temp_cs_name in cs_map:
                cs_map[temp_cs_name].delete()
        except Exception:
            pass

        # Recreate permanent PhaseCentreCS in unit cell design
        permanent_cs_name = "PhaseCentreCS"
        cs_map = {cs.name: cs for cs in source_modeler.coordinate_systems}
        if permanent_cs_name in cs_map:
            try:
                cs_map[permanent_cs_name].delete()
            except Exception:
                pass

        print(f"  Saving permanent '{permanent_cs_name}' coordinate system to unit cell...")
        source_modeler.create_coordinate_system(
            origin=[f"{opt_x_rounded}mm", f"{opt_y_rounded}mm", f"{opt_z_rounded}mm"],
            reference_cs="Global",
            name=permanent_cs_name,
        )

        # Save project to disk
        self.source_design_app.save_project()
        print("  Unit cell design saved.")

        return [opt_x_rounded, opt_y_rounded, opt_z_rounded]

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


if __name__ == "__main__":
    import argparse
    import os

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
        "--use-existing-cs",
        "--use-existing-phase-centre",
        action="store_true",
        default=False,
        help="Use existing PhaseCentreCS in the unit cell without prompting or running Optimetrics",
    )
    parser.add_argument(
        "--centre-freq", "--center-freq", type=float, default=79.0, help="Centre frequency in GHz (default: 79.0)"
    )
    parser.add_argument("--bandwidth", type=float, default=4.0, help="Sweep bandwidth in GHz (default: 4.0)")

    args = parser.parse_args()

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
        non_graphical=True,
    )

    try:
        builder.connect_desktop()
        if args.layout_path:
            print(f"Loading custom layout from: {args.layout_path}")
            elements_list = load_layout_file(args.layout_path)
        else:
            print(f"No custom layout provided. Using default coplanar layout at {args.centre_freq} GHz...")
            elements_list = builder.calculate_default_coplanar_layout(
                transmitter_count=4,
                receiver_count=4,
                operating_frequency_ghz=args.centre_freq,
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
            operating_frequency_ghz=args.centre_freq,
            use_existing_cs=True if (args.use_existing_cs or args.yes) else None,
            overwrite=True if (args.overwrite or args.yes) else None,
            setup_results=setup_results,
            run_simulation=run_simulation,
        )
    finally:
        builder.close()
