#!/usr/bin/env python3
"""
Module 1: HFSS Full-Wave Array Synthesis.
Programmatically builds a planar MIMO antenna array layout on a single PCB substrate in HFSS
using a template-based, boolean-driven assembly workflow.
"""

import os
import sys
import time
import math
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
    ):
        self.project_path = os.path.abspath(project_path)
        self.source_design_name = source_design_name
        self.target_design_name = target_design_name if target_design_name else f"{source_design_name}MimoArray"
        self.pcb_margin_mm = pcb_margin_mm
        self.grpc_port = grpc_port
        self.non_graphical = non_graphical
        
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
            elements.append({
                "label": f"Rx_{rx_idx + 1}",
                "pos": [x_pos, 0.0, 0.0],
                "role": "Rx",
                "polarization": "v",
                "yaw": 0.0,
            })

        # Tx elements along X-axis at y=subarray_spacing
        tx_offset_x = (transmitter_count - 1) * transmitter_spacing_mm / 2.0
        for tx_idx in range(transmitter_count):
            x_pos = tx_idx * transmitter_spacing_mm - tx_offset_x
            elements.append({
                "label": f"Tx_{tx_idx + 1}",
                "pos": [x_pos, subarray_spacing_mm, 0.0],
                "role": "Tx",
                "polarization": "v",
                "yaw": 180.0,
            })

        print(f"Generated coplanar layout at {operating_frequency_ghz} GHz:")
        print(f"  Rx Elements: {receiver_count} (spacing={receiver_spacing_mm:.2f} mm)")
        print(f"  Tx Elements: {transmitter_count} (spacing={transmitter_spacing_mm:.2f} mm)")
        return elements

    def synthesize_array_in_hfss(self, elements, operating_frequency_ghz=79.0):
        """
        Synthesizes the MIMO array in HFSS using a naming-convention boolean assembly:
        1. Classifies objects in the unit cell (Global Layers, Port Sheets, Active Copper, Dummy Cutouts).
        2. Sizes and draws the global continuous substrate and ground layers.
        3. Copies the source templates into the target design exactly ONCE.
        4. Replicates active structures, port sheets, and dummy solids to all coordinates using local duplicate command.
        5. Applies boolean operations using the dummy solids (e.g. Subtract_L2_Ground).
        6. Assigns lumped port excitations to the replicated port sheets.
        """
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
        port_sheets = []    # list of port sheet names
        dummy_solids = {}   # dict of (operation, target_solid) -> list of dummy source names
        active_elements = [] # list of active trace/copper names
        
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

        # --- Check for PhaseCentreCS (British spelling) ---
        offset = [0.0, 0.0, 0.0]
        cs_map = {cs.name: cs for cs in source_modeler.coordinate_systems}
        
        if "PhaseCentreCS" in cs_map:
            cs_obj = cs_map["PhaseCentreCS"]
            # Extract [dx, dy, dz] origin offset of phase centre
            offset = [float(val) for val in cs_obj.origin]
            print(f"\n[INFO] PhaseCentreCS found in unit-cell design. Origin offset: {offset} mm.")
        else:
            print("\n" + "="*80)
            print("[WARNING] PhaseCentreCS was NOT found in the unit-cell design!")
            print("Proceeding without phase centre offset (zero offset) might cause geometric misalignment")
            print("in phase centre spatial diversity calculations.")
            print("="*80 + "\n")
            
            user_input = input("Do you want to proceed without PhaseCentreCS offset? (y/n): ")
            if user_input.strip().lower() not in ["y", "yes", ""]:
                print("[INFO] Aborted by user.")
                self.close()
                sys.exit(0)
            print("[INFO] Proceeding with zero offset [0.0, 0.0, 0.0] mm.")

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
            # This crops out any empty substrate margin and forces the board to terminate directly where the feedlines end.
            y_min_active_uc = None
            y_max_active_uc = None
            for active_name in (active_elements + port_sheets):
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
                    origin=[origin_x, origin_y, z_min],
                    sizes=[board_width, board_height, thickness],
                    name=layer_name,
                    material=material,
                )
                sub_board.transparency = 0.5
                
            elif layer_name.endswith("_Ground"):
                ground_layer_names.append(layer_name)
                # If ground is modeled as infinitely thin sheet (thickness = 0)
                if thickness == 0.0:
                    print(f"  Creating ground plane sheet '{layer_name}' at Z={z_min}")
                    target_modeler.create_rectangle(
                        orientation="XY",
                        origin=[origin_x, origin_y, z_min],
                        sizes=[board_width, board_height],
                        name=layer_name,
                    )
                else:
                    print(f"  Creating ground plane block '{layer_name}' (thickness={thickness:.2f} mm)")
                    target_modeler.create_box(
                        origin=[origin_x, origin_y, z_min],
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

        # Clean up any custom coordinate systems carried over from the unit-cell templates
        cs_names = [cs.name for cs in target_modeler.coordinate_systems]
        if cs_names:
            print(f"Clearing {len(cs_names)} lingering coordinate systems from target design...")
            try:
                target_modeler.delete(cs_names)
            except Exception as e:
                print(f"  Warning: Failed to delete coordinate systems ({e})")

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
            suffix = port_sheet.replace('PortSheet', '')
            if suffix.startswith('_'):
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

        # --- STEP 7: Create Radiation Airbox ---
        print("\nCreating radiation boundary Airbox...")
        # Get overall board bounds in Z direction
        all_z_coords = []
        for layer_info in global_layers.values():
            all_z_coords.extend([layer_info["z_min"], layer_info["z_max"]])
        overall_z_min = min(all_z_coords) if all_z_coords else 0.0
        overall_z_max = max(all_z_coords) if all_z_coords else 1.0

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

        self.target_design_app.save_project()
        print(f"\nMIMO array synthesis successfully completed in design '{self.target_design_name}'.")

    def run_solve(self):
        """Triggers the HFSS simulation setup and solve."""
        if not self.target_design_app:
            raise RuntimeError("Array design not synthesized yet.")
        
        print("Initializing analysis setup...")
        operating_freq = "28GHz"
        if self.source_design_app.setups:
            operating_freq = self.source_design_app.setups[0].props.get("Frequency", "28GHz")
            
        setup = self.target_design_app.create_setup(setup_name="ArraySetup")
        setup.props["Frequency"] = operating_freq
        
        print(f"Solving full-wave HFSS array design at {operating_freq}...")
        self.target_design_app.analyze(setup_name="ArraySetup")

    def export_coupling_s_parameters(self, output_touchstone_path):
        """Exports solved S-parameter coupling matrix to a Touchstone (.sNp) file."""
        if not self.target_design_app:
            raise RuntimeError("Array design is not resolved.")
            
        output_touchstone_path = os.path.abspath(output_touchstone_path)
        print(f"Exporting coupling S-parameters to: {output_touchstone_path}")
        
        self.target_design_app.export_touchstone(
            setup_name="ArraySetup",
            sweep_name="LastSweep",
            filename=output_touchstone_path
        )

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
    if len(sys.argv) < 3:
        print("Usage: python3 hfss_array_builder.py <project_path> <source_design_name>")
        sys.exit(1)

    builder = MimoHfssBuilder(
        project_path=sys.argv[1],
        source_design_name=sys.argv[2],
        non_graphical=True
    )
    try:
        builder.connect_desktop()
        elements_list = builder.calculate_default_coplanar_layout(
            transmitter_count=4,
            receiver_count=4,
            operating_frequency_ghz=79.0,
            subarray_spacing_mm=10.0,
        )
        builder.synthesize_array_in_hfss(elements_list, operating_frequency_ghz=79.0)
    finally:
        builder.close()
