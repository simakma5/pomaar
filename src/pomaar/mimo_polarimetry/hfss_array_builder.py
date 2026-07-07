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

    def synthesize_array_in_hfss(self, elements, operating_frequency_ghz=None):
        """
        Synthesizes the MIMO array in HFSS using a naming-convention boolean assembly:
        1. Classifies objects in the unit cell (Global Layers, Port Sheets, Active Copper, Dummy Cutouts).
        2. Sizes and draws the global continuous substrate and ground layers.
        3. Copies the source templates into the target design exactly ONCE.
        4. Replicates active structures, port sheets, and dummy solids to all coordinates using local duplicate command.
        5. Applies boolean operations using the dummy solids (e.g. Subtract_L2_Ground).
        6. Assigns lumped port excitations to the replicated port sheets.
        """
        if operating_frequency_ghz is not None:
            self.centre_frequency_ghz = operating_frequency_ghz
        else:
            operating_frequency_ghz = self.centre_frequency_ghz

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

        # --- STEP 1: Classify Unit-Cell Objects ---
        # Perform object scan while source design is active
        all_object_names = source_modeler.object_names

        global_layers = {}  # dict of layer_name -> {material, z_min, z_max}
        port_sheets = []  # list of port sheet names
        dummy_solids = {}  # dict of (operation, target_solid) -> list of dummy source names
        active_elements = []  # list of active trace/copper names

        for obj_name in all_object_names:
            # 1. First, check for dummy boolean solids (e.g. Subtract_L12_Ground)
            if "_" in obj_name and obj_name.split("_")[0] in ["Subtract", "Unite"]:
                tokens = obj_name.split("_")
                operation = tokens[0]
                target_solid = "_".join(tokens[1:])
                key = (operation, target_solid)
                if key not in dummy_solids:
                    dummy_solids[key] = []
                dummy_solids[key].append(obj_name)

            # 2. Second, check for global layers (Substrate or Ground)
            elif obj_name.endswith("_Substrate") or obj_name.endswith("_Ground"):
                source_obj = source_modeler.get_object_from_name(obj_name)
                bbox = source_obj.bounding_box
                z_min, z_max = float(bbox[2]), float(bbox[5])
                global_layers[obj_name] = {
                    "material": source_obj.material_name,
                    "z_min": z_min,
                    "z_max": z_max,
                }

            # 3. Third, check for port sheets (PortSheet1, PortSheet2, etc.)
            elif obj_name.startswith("PortSheet"):
                port_sheets.append(obj_name)

            # 4. Fourth, regular active elements (Patches, Feeds, etc.)
            else:
                source_obj = source_modeler.get_object_from_name(obj_name)
                # Skip vacuum objects (like BoundingBox, RadiatingSurface) to avoid clutter
                if source_obj and source_obj.material_name.lower() == "vacuum":
                    print(f"  Skipping vacuum object: '{obj_name}'")
                    continue
                active_elements.append(obj_name)

        print("\nUnit-cell structure classified:")
        print(f"  Global Layers:    {list(global_layers.keys())}")
        print(f"  Port Sheets:      {port_sheets}")
        print(f"  Dummy Cutouts:    {list(dummy_solids.keys())}")
        print(f"  Active Elements:  {active_elements}")

        # --- Check/Compute PhaseCentreCS ---
        offset = [0.0, 0.0, 0.0]
        cs_map = {cs.name: cs for cs in source_modeler.coordinate_systems}

        run_opt = False
        if "PhaseCentreCS" in cs_map:
            cs_obj = cs_map["PhaseCentreCS"]
            try:
                # Direct evaluation to bypass PyAEDT temp_var post-processing assignment bug
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
            if sys.stdin.isatty():
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
            if sys.stdin.isatty():
                user_input = input("Do you want to run Optimetrics to calculate the Phase Centre? (y/n) [default: y]: ")
            else:
                user_input = "y"
            if user_input.strip().lower() not in ["n", "no"]:
                run_opt = True
            else:
                if sys.stdin.isatty():
                    user_input2 = input("Proceed without PhaseCentreCS offset? (y/n) [default: y]: ")
                else:
                    user_input2 = "y"
                if user_input2.strip().lower() in ["n", "no"]:
                    print("[INFO] Aborted by user.")
                    self.close()
                    sys.exit(0)
                print("[INFO] Proceeding with zero offset [0.0, 0.0, 0.0] mm.")

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
            if sys.stdin.isatty():
                user_input = input(f"Design '{self.target_design_name}' already exists. Do you want to clear it and overwrite? (y/n) [default: y]: ")
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

            # Clear existing 3D bodies
            target_modeler = self.target_design_app.modeler
            all_objs = list(target_modeler.object_names)
            if all_objs:
                print(f"  Clearing {len(all_objs)} existing objects...")
                target_modeler.delete(all_objs)

            # Clear coordinate systems
            cs_names = [cs.name for cs in target_modeler.coordinate_systems]
            if cs_names:
                print(f"  Clearing {len(cs_names)} coordinate systems...")
                target_modeler.delete(cs_names)

            # Setup and sweep purging will be handled at the end of synthesis in configure_simulation_setup
        else:
            print(f"Creating new HFSS array design '{self.target_design_name}'...")
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

        # --- STEP 3: Size & Draw Global Contiguous Board ---
        # Sort layers by prefix (L1, L12, L2...) to ensure they are created in stackup order
        sorted_layers = sorted(global_layers.keys(), key=lambda name: name.split("_")[0])
        ground_layer_names = []

        # Track overall board bounds for Airbox radiation boundary creation
        overall_x_min = None
        overall_x_max = None
        overall_y_min = None
        overall_y_max = None

        print("\nConstructing global contiguous board layers:")
        for layer_name in sorted_layers:
            layer_info = global_layers[layer_name]
            z_min = layer_info["z_min"]
            z_max = layer_info["z_max"]
            thickness = z_max - z_min
            material = layer_info["material"]

            # Query the unit-cell object to get its exact bounding box (for X-bounds)
            layer_obj = self.source_design_app.modeler.get_object_from_name(layer_name)
            bbox = layer_obj.bounding_box
            x_min_uc, x_max_uc = float(bbox[0]), float(bbox[3])

            # Query all active elements and port sheets to get the Y-bounds (longitudinal bounds)
            # This crops out empty substrate margin and forces the board to terminate directly where the feedlines end.
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

            # Fallback to substrate bounds if no active elements are present
            if y_min_active_uc is None:
                y_min_active_uc = float(bbox[1])
            if y_max_active_uc is None:
                y_max_active_uc = float(bbox[4])

            # Define the 4 corners: X matches the substrate edges, Y matches the feedline edges
            corners = [
                (x_min_uc, y_min_active_uc),
                (x_max_uc, y_min_active_uc),
                (x_min_uc, y_max_active_uc),
                (x_max_uc, y_max_active_uc),
            ]

            # Calculate global bounding box by union of all rotated/translated footprints
            all_x_glob = []
            all_y_glob = []
            for element in elements:
                pos = element["pos"]
                yaw_deg = element.get("yaw", 0.0)
                rad = math.radians(yaw_deg)
                cos_val = math.cos(rad)
                sin_val = math.sin(rad)

                for cx, cy in corners:
                    # Rotate corner around origin
                    rx = cx * cos_val - cy * sin_val
                    ry = cx * sin_val + cy * cos_val
                    # Translate to element position
                    all_x_glob.append(rx + pos[0])
                    all_y_glob.append(ry + pos[1])

            global_x_min = min(all_x_glob) - self.pcb_margin_mm
            global_x_max = max(all_x_glob) + self.pcb_margin_mm
            global_y_min = min(all_y_glob) - self.pcb_margin_mm
            global_y_max = max(all_y_glob) + self.pcb_margin_mm

            # Track overall bounds across all board layers
            if overall_x_min is None or global_x_min < overall_x_min:
                overall_x_min = global_x_min
            if overall_x_max is None or global_x_max > overall_x_max:
                overall_x_max = global_x_max
            if overall_y_min is None or global_y_min < overall_y_min:
                overall_y_min = global_y_min
            if overall_y_max is None or global_y_max > overall_y_max:
                overall_y_max = global_y_max

            board_width = global_x_max - global_x_min
            board_height = global_y_max - global_y_min
            origin_x = global_x_min
            origin_y = global_y_min

            if layer_name.endswith("_Substrate"):
                print(f"  Creating substrate '{layer_name}' (material={material}, thickness={thickness:.2f} mm)")
                sub_board = target_modeler.create_box(
                    origin=[origin_x, origin_y, z_min - offset[2]],
                    sizes=[board_width, board_height, thickness],
                    name=layer_name,
                    material=material,
                )
                sub_board.transparency = 0.5

            elif layer_name.endswith("_Ground"):
                ground_layer_names.append(layer_name)
                # If ground is modeled as infinitely thin sheet (thickness = 0)
                if thickness == 0.0:
                    print(f"  Creating ground plane sheet '{layer_name}' at Z={z_min - offset[2]}")
                    target_modeler.create_rectangle(
                        orientation="XY",
                        origin=[origin_x, origin_y, z_min - offset[2]],
                        sizes=[board_width, board_height],
                        name=layer_name,
                    )
                else:
                    print(f"  Creating ground plane block '{layer_name}' (thickness={thickness:.2f} mm)")
                    target_modeler.create_box(
                        origin=[origin_x, origin_y, z_min - offset[2]],
                        sizes=[board_width, board_height, thickness],
                        name=layer_name,
                        material=material,
                    )

        # --- STEP 4: Replicate and Translate Elements ---
        all_replicate_sources = active_elements + port_sheets
        for dummy_list in dummy_solids.values():
            all_replicate_sources.extend(dummy_list)

        replicated_dummy_mapping = {key: [] for key in dummy_solids.keys()}
        replicated_ports_list = []

        # Optimization: Copy the templates from the source design to the target design exactly ONCE.
        print("\nCopying template geometries to target design...")
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

        print("Replicating structures to array grid...")
        for element in elements:
            label = element["label"]
            pos = element["pos"]
            polarization = element.get("polarization", "v").lower()
            element_yaw = element.get("yaw", 0.0)

            polarization_variants = []
            if polarization.lower() == "v" or polarization.lower() == "both":
                polarization_variants.append(("_V", 0.0))  # Vertical (no rotation)
            if polarization.lower() == "h" or polarization.lower() == "both":
                polarization_variants.append(("_H", 90.0))  # Horizontal (90-deg Z rotation)

            for variant_suffix, rotation_yaw in polarization_variants:
                variant_label = f"{label}{variant_suffix}"
                print(f"  Replicating element variant: {variant_label} to position {pos} mm...")

                # Duplicate all templates locally using duplicate_along_line with a dummy Z offset of 1.0 mm.
                # This completely bypasses the X11 clipboard copy/paste mechanism to avoid hangs in headless containers.
                success, pasted_names = target_modeler.duplicate_along_line(
                    assignment=all_replicate_sources,
                    vector=[0, 0, 1.0],
                    clones=2,
                )
                if not success:
                    raise RuntimeError(f"Failed to duplicate template objects for variant {variant_label}")

                # Move back to origin (offset the Z translation)
                target_modeler.move(
                    assignment=pasted_names,
                    vector=[0, 0, -1.0],
                )

                # Rename duplicated objects and track ports/dummy solids
                renamed_objs = []
                for pasted_name in pasted_names:
                    # Strip numerical suffixes added by AEDT duplicate if any (e.g. L1_Patch_1 -> L1_Patch)
                    base_name = pasted_name
                    for source_name in all_replicate_sources:
                        if pasted_name.startswith(source_name):
                            base_name = source_name
                            break

                    new_name = f"{base_name}_{variant_label}"
                    target_modeler.get_object_from_name(pasted_name).name = new_name
                    renamed_objs.append(new_name)

                    # Track port sheets
                    if base_name in port_sheets:
                        replicated_ports_list.append(new_name)

                    # Track dummy solids
                    for key, dummy_src_list in dummy_solids.items():
                        if base_name in dummy_src_list:
                            replicated_dummy_mapping[key].append(new_name)

                # Rotate variant if needed (centred at the origin)
                total_rotation = element_yaw + rotation_yaw
                if total_rotation != 0.0:
                    target_modeler.rotate(
                        assignment=renamed_objs,
                        axis=Axis.Z,
                        angle=total_rotation,
                    )

                # Calculate the rotated phase centre offset
                rad = math.radians(total_rotation)
                cos_val = math.cos(rad)
                sin_val = math.sin(rad)

                # Rotate offset around Z axis
                rot_dx = offset[0] * cos_val - offset[1] * sin_val
                rot_dy = offset[0] * sin_val + offset[1] * cos_val
                rot_dz = offset[2]

                # Move variant to final position (compensating for the phase centre offset)
                target_modeler.move(
                    assignment=renamed_objs,
                    vector=[pos[0] - rot_dx, pos[1] - rot_dy, pos[2] - rot_dz],
                )

        # Clean up the original template geometries at the origin of the target design
        print("Cleaning up template geometries...")
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

        # --- STEP 6: Assign Lumped Port Excitations ---
        print("\nAssigning port excitations...")
        # Use first identified ground plane as the default reference ground
        default_ground = ground_layer_names[0] if ground_layer_names else "Ground_Plane"

        for port_sheet in replicated_ports_list:
            suffix = port_sheet.replace("PortSheet", "")
            if suffix.startswith("_"):
                suffix = suffix[1:]
            port_name = f"Port_{suffix}"

            print(f"  Assigning lumped port to sheet '{port_sheet}' referencing '{default_ground}'")
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

        # Calculate clearance using lambda_0 / 4 (quarter wavelength) rule of thumb
        speed_of_light_mm_s = 2.99792458e11
        wavelength_mm = speed_of_light_mm_s / (operating_frequency_ghz * 1e9)
        airbox_clearance_mm = 0.25 * wavelength_mm

        airbox_x_min = overall_x_min - airbox_clearance_mm
        airbox_x_max = overall_x_max + airbox_clearance_mm
        airbox_y_min = overall_y_min - airbox_clearance_mm
        airbox_y_max = overall_y_max + airbox_clearance_mm
        airbox_z_min = overall_z_min - airbox_clearance_mm
        airbox_z_max = overall_z_max + airbox_clearance_mm

        airbox_width = airbox_x_max - airbox_x_min
        airbox_height = airbox_y_max - airbox_y_min
        airbox_thickness = airbox_z_max - airbox_z_min

        print(f"  Creating Airbox solid (clearance={airbox_clearance_mm:.2f} mm)")
        airbox_obj = target_modeler.create_box(
            origin=[airbox_x_min, airbox_y_min, airbox_z_min],
            sizes=[airbox_width, airbox_height, airbox_thickness],
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

        print("  Hiding Airbox in layout view (making it fully transparent)...")
        try:
            airbox_obj.transparency = 1.0
        except Exception as e:
            print(f"  Warning: Failed to set Airbox transparency ({e})")

        # Configure simulation setup and sweep in the target design
        self.configure_simulation_setup()

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
            f"Configuring interpolating frequency sweep '{sweep_name}' ({start_freq:.2f} GHz - {end_freq:.2f} GHz, 401 points)..."
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

    def run_solve(self):
        """Triggers the HFSS simulation setup and solve."""
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")

        # Ensure the setup is configured
        self.configure_simulation_setup()

        print("Solving full-wave HFSS array design using 'ArraySetup'...")
        self.target_design_app.analyze(setup_name="ArraySetup")

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
        opt_setup = self.source_design_app.optimizations.add(
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


if __name__ == "__main__":
    import argparse
    import json
    import os

    parser = argparse.ArgumentParser(description="HFSS Full-Wave Array Synthesis CLI.")
    parser.add_argument("project_path", help="Path to the AEDT project file (.aedt)")
    parser.add_argument("source_design_name", help="Name of the unit-cell source design")
    parser.add_argument(
        "layout_json_path", nargs="?", default=None, help="Optional path to custom elements layout JSON file"
    )
    parser.add_argument(
        "--centre-freq", "--center-freq", type=float, default=79.0, help="Centre frequency in GHz (default: 79.0)"
    )
    parser.add_argument("--bandwidth", type=float, default=4.0, help="Sweep bandwidth in GHz (default: 4.0)")

    args = parser.parse_args()

    # Determine target design name from layout file if provided
    target_design_name = None
    if args.layout_json_path:
        base_name = os.path.splitext(os.path.basename(args.layout_json_path))[0]
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
        if args.layout_json_path:
            print(f"Loading custom layout from: {args.layout_json_path}")
            with open(args.layout_json_path, "r") as f:
                elements_list = json.load(f)
        else:
            print(f"No custom layout provided. Using default coplanar layout at {args.centre_freq} GHz...")
            elements_list = builder.calculate_default_coplanar_layout(
                transmitter_count=4,
                receiver_count=4,
                operating_frequency_ghz=args.centre_freq,
                subarray_spacing_mm=10.0,
            )
        builder.synthesize_array_in_hfss(elements_list, operating_frequency_ghz=args.centre_freq)
    finally:
        builder.close()
