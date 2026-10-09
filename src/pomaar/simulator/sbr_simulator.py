#!/usr/bin/env python3
"""
Module 2: SBR+ Target Solver.
Builds an SBR+ design around a canonical target (sphere, dihedral or trihedral) and links a
solved HFSS antenna design from the same project as its source, then sets up the Tx/Rx pairs,
the analysis setup and quick-look S-parameter reports.
"""

import argparse
import math
import os
import re
import sys

from ansys.aedt.core import Desktop, Hfss
from ansys.aedt.core.generic.constants import Plane
from ansys.aedt.core.generic.general_methods import grpc_active_sessions

from pomaar.simulator.hfss_array_builder import (
    DEFAULT_CENTRE_FREQUENCY_GHZ,
    SPEED_OF_LIGHT_MM_S,
    BuildAbortedError,
    NativeAedtRefusedError,
    _confirm,
    _length_mm,
    _parse_frequency_ghz,
)

# 'none' builds the same scene without a target; its solve is the background (direct Tx -> Rx
# coupling between the proxy ports) to subtract coherently from target runs
TARGET_TYPES = ("sphere", "dihedral", "trihedral", "none")
FIELD_TYPES = ("farfield", "nearfield")

# Defaults for 79 GHz (lambda = 3.8 mm). 'size' is the sphere radius or the dihedral/trihedral
# edge length. A 40 mm target has 2 D^2 / lambda <= 3.4 m for all three types (D = 2a, a*sqrt(3)
# and a*sqrt(2)), so a 4 m range keeps it in its own far field.
DEFAULT_TARGET_DISTANCE = "4meter"
DEFAULT_TARGET_SIZE = "40mm"

ANTENNA_CS = "AntennaCS"
TARGET_CS = "TargetCS"
SETUP_NAME = "Setup"
SWEEP_NAME = "Sweep"

# Boundaries whose objects only bound the full-wave solution space; they are left out of the
# linked antenna's Model Visualization. The builder's radiation box is called 'Airbox'.
RADIATION_BOUNDARY_TYPES = ("Radiation", "FE-BI", "PML", "Hybrid")
RADIATION_OBJECT_NAMES = ("Airbox",)

# SBR+ splits the solution 'Setup : Sweep' into '_Incident' (the direct Tx -> Rx coupling between
# the proxy ports) and '_Scattered' (the target return); target signatures use the latter
SCATTERED_SWEEP_NAME = f"{SWEEP_NAME}_Scattered"

REPORT_S_DB = "S-parameters dB"
REPORT_S_DB_SCATTERED = "S-parameters dB (scattered)"
REPORT_COPOL_IMBALANCE = "Co-pol VV/HH"
REPORT_CROSSPOL_ISOLATION = "Cross-pol isolation"


def _parse_port(port_name):
    """
    Classifies a source port name as (role, polarization, element key), e.g. 'TxH' ->
    ('Tx', 'H', 'Tx') and 'Port_Tx_1_H' -> ('Tx', 'H', 'Port_Tx_1_'); None if it is neither.
    The element key pairs the H and V ports of one dual-polarised element.
    """
    role = re.search(r"(Tx|Rx)", port_name)
    pol = re.search(r"([HV])$", port_name)
    if not role or not pol:
        return None
    return role.group(1), pol.group(1), port_name[:-1]


def _ask_field_type(preset=None):
    """Resolves the linked antenna field type: a preset wins, then a prompt, then far field."""
    if preset is not None:
        if preset not in FIELD_TYPES:
            raise ValueError(f"field_type must be one of {FIELD_TYPES}, not '{preset}'")
        return preset
    if not sys.stdin.isatty():
        return "farfield"
    answer = (
        input(
            "Link the antenna with far-field results (f, standard) or near-field results "
            "(n, e.g. bumper/radome integration)? [default: f]: "
        )
        .strip()
        .lower()
    )
    return "nearfield" if answer in ["n", "near", "nearfield"] else "farfield"


class SbrSimulationManager:
    """
    Builds an SBR+ design that links a solved HFSS antenna design (same project) as the source and
    places a canonical PEC target on a sphere around the antenna.

    Geometry and variables of the SBR+ design:
    * AntennaCS sits at the global origin; the linked antenna is placed in it (boresight +z).
    * TargetCS is defined in AntennaCS at the spherical position
      (targetDistance, targetAzimuth, targetElevation), with targetElevation the polar angle off
      boresight: X = d cos(az) sin(el), Y = d sin(az) sin(el), Z = d cos(el).
    * TargetCS is rotated by Euler angles ZYZ (targetAzimuth, targetElevation,
      targetRoll - targetAzimuth), which turn the local z axis radially outward without twisting
      it: the target always faces the antenna and targetRoll is a pure rotation about the line of
      sight. Edit the TargetCS angles directly for other orientations.
    * All targets face the antenna along local -z, scale with targetSize and carry a Perfect E
      boundary: the PEC sphere (radius) is centred at TargetCS; the dihedral seam (seam length =
      sheet depth = targetSize, along local x, i.e. parallel to H) and the trihedral corner (edge
      length targetSize) of the sheet targets sit at TargetCS.
    """

    def __init__(
        self,
        project_path,
        source_design_name,
        sbr_design_name=None,
        grpc_port=None,
        non_graphical=True,
    ):
        self.project_path = os.path.abspath(project_path)
        self.source_design_name = source_design_name
        self.sbr_design_name = sbr_design_name if sbr_design_name else f"Linked{source_design_name}"
        self.grpc_port = grpc_port
        self.non_graphical = non_graphical

        self.desktop_session = None
        self.source_app = None
        self.sbr_app = None
        self.is_new_desktop = False

        self.antenna_instance_name = None
        self.source_solution = None
        self.tx_ports = []
        self.rx_ports = []
        self.target_object_names = []

    def connect_desktop(self):
        """Starts or connects to an existing AEDT session via gRPC (see MimoHfssBuilder)."""
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
            print(f"Failed to connect ({err}). Starting new session...")
            self.desktop_session = Desktop(
                port=self.grpc_port,
                new_desktop=True,
                non_graphical=self.non_graphical,
                close_on_exit=True,
            )
            self.is_new_desktop = True

    # ------------------------------------------------------------------------------------------
    # Full flow
    # ------------------------------------------------------------------------------------------

    def setup_sbr_design(
        self,
        target_type="dihedral",
        target_size=DEFAULT_TARGET_SIZE,
        target_roughness=0.0,
        target_height_deviation=0.0,
        target_distance=DEFAULT_TARGET_DISTANCE,
        target_azimuth="0deg",
        target_elevation="0deg",
        target_roll="0deg",
        field_type=None,
        mapped_variables=None,
        source_solution=None,
        ray_density=4,
        max_bounces=5,
        overwrite=None,
        simulate=None,
    ):
        """
        Builds the complete SBR+ design and optionally solves it.

        Parameters
        ----------
        target_type : str
            'sphere', 'dihedral' (default) or 'trihedral'.
        target_size : str or float
            Sphere radius or dihedral/trihedral edge length (floats are mm).
        target_roughness, target_height_deviation : float, optional
            SBR+ surface roughness and surface height standard deviation (mm) of the target's
            Perfect E boundary; both 0 (smooth) by default.
        target_distance, target_azimuth, target_elevation, target_roll : str or float
            Initial values of the TargetCS design variables (floats are m and deg).
        field_type : str, optional
            'farfield' or 'nearfield'; asked interactively when None.
        mapped_variables : list of str, optional
            Independent variables of the HFSS design that are mapped to SBR+ design variables
            of the same name (created from the HFSS value when missing).
        source_solution : str, optional
            Linked HFSS solution; defaults to the first frequency sweep of the HFSS design.
        overwrite, simulate : bool, optional
            Presets for the 'replace existing design' and 'start simulation' prompts.
        """
        if not self.desktop_session:
            self.connect_desktop()

        self.create_sbr_design(overwrite=overwrite)
        self.create_coordinate_systems(
            target_distance=target_distance,
            target_azimuth=target_azimuth,
            target_elevation=target_elevation,
            target_roll=target_roll,
        )
        self.create_target(
            target_type=target_type,
            size=target_size,
            roughness=target_roughness,
            height_deviation=target_height_deviation,
        )
        self.link_antenna(
            field_type=field_type,
            mapped_variables=mapped_variables,
            source_solution=source_solution,
        )
        self.assign_tx_rx()
        self.create_analysis_setup(ray_density=ray_density, max_bounces=max_bounces)
        self.check_target_far_field(target_type)
        self.create_reports()
        self.sbr_app.save_project()
        self.run_solve(simulate=simulate)

    # ------------------------------------------------------------------------------------------
    # Steps
    # ------------------------------------------------------------------------------------------

    def create_sbr_design(self, overwrite=None):
        """Opens the HFSS source design and (re)creates the SBR+ design next to it."""
        print(f"Loading project: {self.project_path}")
        self.source_app = Hfss(
            project=self.project_path,
            design=self.source_design_name,
            new_desktop=False,
            close_on_exit=False,
        )
        if self.source_app.solution_type == "SBR+":
            raise ValueError(f"Source design '{self.source_design_name}' is itself an SBR+ design.")

        if self.sbr_design_name in self.source_app.design_list:
            question = (
                f"SBR+ design '{self.sbr_design_name}' already exists. Delete and rebuild it?"
            )
            if not _confirm(question, preset=overwrite, default=False):
                raise BuildAbortedError(f"Kept existing SBR+ design '{self.sbr_design_name}'.")
            # AEDT refuses to reuse a deleted design's name until the project is saved, so the old
            # design is renamed out of the way before it is deleted
            print(f"Deleting existing SBR+ design '{self.sbr_design_name}'...")
            oproject = self.source_app.oproject
            retired_name = f"{self.sbr_design_name}Deleted"
            oproject.GetDesign(self.sbr_design_name).RenameDesignInstance(
                self.sbr_design_name, retired_name
            )
            oproject.DeleteDesign(retired_name)

        # Native calls: PyAEDT's insert_design/delete_design re-point source_app to another design
        print(f"Creating SBR+ design '{self.sbr_design_name}'...")
        self.source_app.oproject.InsertDesign("HFSS", self.sbr_design_name, "SBR+", "")
        if self.sbr_design_name not in self.source_app.oproject.GetTopDesignList():
            raise RuntimeError(f"AEDT did not create SBR+ design '{self.sbr_design_name}'.")
        self.sbr_app = Hfss(
            project=self.source_app.project_name,
            design=self.sbr_design_name,
            new_desktop=False,
            close_on_exit=False,
        )
        self.sbr_app.modeler.model_units = "mm"

    def create_coordinate_systems(
        self,
        target_distance=DEFAULT_TARGET_DISTANCE,
        target_azimuth="0deg",
        target_elevation="0deg",
        target_roll="0deg",
    ):
        """(a) AntennaCS at the global origin; (b) spherically placed, ZYZ-rotated TargetCS."""
        app = self.sbr_app
        variables = [
            ("targetDistance", _with_unit(target_distance, "meter")),
            ("targetAzimuth", _with_unit(target_azimuth, "deg")),
            ("targetElevation", _with_unit(target_elevation, "deg")),
            ("targetRoll", _with_unit(target_roll, "deg")),
        ]
        for name, value in variables:
            app[name] = value

        editor = app.modeler.oeditor
        print(f"Creating {ANTENNA_CS} at the global origin...")
        app.modeler.set_working_coordinate_system("Global")
        editor.CreateRelativeCS(
            _euler_zyz_cs_args(["0mm", "0mm", "0mm"], ["0deg", "0deg", "0deg"]),
            ["NAME:Attributes", "Name:=", ANTENNA_CS],
        )

        # CreateRelativeCS references the working CS, so TargetCS is created from AntennaCS
        print(f"Creating {TARGET_CS} in {ANTENNA_CS} (spherical position, Euler ZYZ rotation)...")
        app.modeler.set_working_coordinate_system(ANTENNA_CS)
        editor.CreateRelativeCS(
            _euler_zyz_cs_args(
                [
                    "targetDistance*cos(targetAzimuth)*sin(targetElevation)",
                    "targetDistance*sin(targetAzimuth)*sin(targetElevation)",
                    "targetDistance*cos(targetElevation)",
                ],
                # Rz(az) Ry(el) Rz(roll - az): local z radially outward, no twist at roll = 0
                ["targetAzimuth", "targetElevation", "targetRoll - targetAzimuth"],
            ),
            ["NAME:Attributes", "Name:=", TARGET_CS],
        )
        app.modeler.set_working_coordinate_system("Global")

    def create_target(
        self, target_type="dihedral", size=DEFAULT_TARGET_SIZE, roughness=0.0, height_deviation=0.0
    ):
        """
        (c) Target in TargetCS, scaled by targetSize, with a Perfect E boundary (SBR+ surface
        roughness and height standard deviation, both 0 by default): a PEC sphere, or a dihedral
        or trihedral of sheets. One sheet is drawn and the others are rotation copies.
        """
        target_type = target_type.lower()
        if target_type not in TARGET_TYPES:
            raise ValueError(f"target_type must be one of {TARGET_TYPES}, not '{target_type}'")
        app = self.sbr_app
        modeler = app.modeler
        app["targetSize"] = _with_unit(size, "mm")
        if target_type == "none":
            print("No target: this design solves the background (direct Tx/Rx coupling) only.")
            self.target_object_names = []
            return

        print(f"Creating {target_type} target (targetSize = {app['targetSize']})...")
        modeler.set_working_coordinate_system(TARGET_CS)
        a = "targetSize"
        if target_type == "sphere":
            target = modeler.create_sphere([0, 0, 0], a, name="TargetSphere", material="pec")
        else:
            if target_type == "dihedral":
                # Square sheet on the local xy plane with its seam edge on the x axis; tilted by
                # -45 deg about x it runs from the seam towards the antenna (-z), and its copy at
                # 180 deg about the symmetry axis z closes the 90 deg fold
                target = modeler.create_rectangle(
                    Plane.XY, [f"-{a}/2", "0mm", "0mm"], [a, a], name="TargetDihedral"
                )
                modeler.rotate(target.name, "X", -45)
                axis, angle, clones = "Z", 180, 2
            else:
                # Triangular sheet spanned by the orthonormal corner edges e1, e2 (e1 + e2 + e3
                # along -z, so the symmetry axis is the local z axis and points at the antenna);
                # its copies at 120 deg about z are the other two faces
                e1 = [f"{a}/sqrt(2)", f"-{a}/sqrt(6)", f"-{a}/sqrt(3)"]
                e2 = [f"-{a}/sqrt(2)", f"-{a}/sqrt(6)", f"-{a}/sqrt(3)"]
                target = modeler.create_polyline(
                    [["0mm", "0mm", "0mm"], e1, e2],
                    close_surface=True,
                    cover_surface=True,
                    name="TargetTrihedral",
                )
                axis, angle, clones = "Z", 120, 3
            success, copies = modeler.duplicate_around_axis(target.name, axis, angle, clones)
            if not success:
                raise RuntimeError(f"Failed to copy the {target_type} plate.")
            modeler.unite([target.name] + list(copies))
        modeler.set_working_coordinate_system("Global")

        # Sheets carry no material; the Perfect E boundary makes every target a PEC surface
        app.assign_perfect_e(
            target.name,
            height_deviation=height_deviation,
            roughness=roughness,
            name=f"{target.name}PerfectE",
        )
        self.target_object_names = [target.name]

    def link_antenna(self, field_type=None, mapped_variables=None, source_solution=None):
        """
        (d) Links the HFSS design (same project) as an antenna in AntennaCS with 'Simulate source
        design as needed', the chosen variable mapping and Model Visualization of all objects
        except the radiation boundary.
        """
        app = self.sbr_app
        src = self.source_app
        field_type = _ask_field_type(field_type)
        self.source_solution = source_solution or self._default_source_solution()

        independent = list(src.variable_manager.independent_design_variable_names)
        expressions = {v: src.variable_manager.variables[v].expression for v in independent}
        mapped = list(mapped_variables or [])
        for name in mapped:
            if name not in expressions:
                raise ValueError(
                    f"Cannot map '{name}': not an independent design variable of "
                    f"'{self.source_design_name}' ({', '.join(independent)})."
                )
            if name not in app.variable_manager.variables:
                app[name] = expressions[name]
        params = ["NAME:Params"]
        for name, expression in expressions.items():
            params += [f"{name}:=", name if name in mapped else expression]

        visualization = self._visualization_objects()
        provider = [
            "NAME:NativeComponentDefinitionProvider",
            "Type:=", "Linked Antenna",
            "Unit:=", app.modeler.model_units,
            "Version:=", 0,
            "Is Parametric Array:=", False,
            "Project:=", "This Project*",
            "Product:=", "HFSS",
            "Design:=", self.source_design_name,
            "Soln:=", self.source_solution,
            params,
            "ForceSourceToSolve:=", True,
            "PreservePartnerSoln:=", False,
            "PathRelativeTo:=", "TargetProject",
            "FieldType:=", field_type,
            "UseCompositePort:=", False,
            ["NAME:SourceBlockageStructure", "NonModelObject:=", []],
            "VisualizationObjects:=", visualization,
        ]  # fmt: skip
        if field_type == "nearfield":
            # PyAEDT's near-field defaults (create_sbr_linked_antenna)
            provider += [
                "UseGlobalCurrentSrcOption:=", True,
                "Current Source Conformance:=", "Disable",
                "Thin Sources:=", True,
                "Power Fraction:=", "0.95",
            ]  # fmt: skip

        component = self.source_design_name
        args = [
            "NAME:InsertNativeComponentData",
            "TargetCS:=", ANTENNA_CS,
            "SubmodelDefinitionName:=", component,
            ["NAME:ComponentPriorityLists"],
            "NextUniqueID:=", 0,
            "MoveBackwards:=", False,
            "DatasetType:=", "ComponentDatasetType",
            ["NAME:DatasetDefinitions"],
            [
                "NAME:BasicComponentInfo",
                "ComponentName:=", component,
                "Company:=", "",
                "Company URL:=", "",
                "Model Number:=", "",
                "Help URL:=", "",
                "Version:=", "1.0",
                "Notes:=", "",
                "IconType:=", "Linked Antenna",
            ],
            ["NAME:GeometryDefinitionParameters", ["NAME:VariableOrders"]],
            ["NAME:DesignDefinitionParameters", ["NAME:VariableOrders"]],
            ["NAME:MaterialDefinitionParameters", ["NAME:VariableOrders"]],
            "DefReferenceCSID:=", 1,
            "MapInstanceParameters:=", "DesignVariable",
            "UniqueDefinitionIdentifier:=", "",
            "OriginFilePath:=", "",
            "IsLocal:=", False,
            "ChecksumString:=", "",
            "ChecksumHistory:=", [],
            "VersionHistory:=", [],
            provider,
            [
                "NAME:InstanceParameters",
                # Each mapped component parameter takes the SBR+ design variable of the same name
                "GeometryParameters:=", " ".join(f"{name}='{name}'" for name in mapped),
                "MaterialParameters:=", "",
                "DesignParameters:=", "",
            ],
        ]  # fmt: skip

        print(
            f"Linking '{self.source_design_name}' ({self.source_solution}, {field_type}) in "
            f"{ANTENNA_CS}; mapped variables: {', '.join(mapped) or 'none'}..."
        )
        self.antenna_instance_name = app.modeler.oeditor.InsertNativeComponent(args)
        print(f"  Linked antenna '{self.antenna_instance_name}'; visualized: {visualization}")

    def assign_tx_rx(self):
        """(e) Every Tx port of the linked antenna transmits to every Rx port (Select Tx/Rx)."""
        tx_ports, rx_ports = [], []
        for port in self.source_app.ports:
            parsed = _parse_port(port)
            if parsed is None:
                print(f"  Warning: cannot tell Tx/Rx and H/V from port '{port}'; it is skipped.")
                continue
            excitation = f"{self.antenna_instance_name}_{port}"
            (tx_ports if parsed[0] == "Tx" else rx_ports).append(excitation)
        if not tx_ports or not rx_ports:
            raise ValueError(f"Linked antenna needs Tx and Rx ports, got {self.source_app.ports}.")

        # 'Rx Antennas' is one comma-separated string; a Python list leaves the matrix empty
        args = ["NAME:SBRTxRxSettings"]
        for index, tx in enumerate(tx_ports):
            args.append(
                [
                    f"NAME:Tx/Rx List {index}",
                    "Tx Antenna:=",
                    tx,
                    "Rx Antennas:=",
                    ",".join(rx_ports),
                ]
            )
        self.sbr_app.oboundary.SetSBRTxRxSettings(args)
        self.tx_ports, self.rx_ports = tx_ports, rx_ports
        print(f"Tx: {tx_ports}\nRx: {rx_ports}")

    def create_analysis_setup(self, ray_density=4, max_bounces=5):
        """(f) SBR+ setup whose frequency sweep copies the linked HFSS sweep."""
        sweep = self._source_sweep()
        sweep_args = [f"NAME:{SWEEP_NAME}"]
        for key, value in sweep.items():
            sweep_args += [f"{key}:=", value]
        print(
            f"Creating SBR+ setup '{SETUP_NAME}' ({sweep['RangeType']} {sweep['RangeStart']} - "
            f"{sweep['RangeEnd']}, ray density {ray_density}, {max_bounces} bounces, "
            "creeping waves)..."
        )
        self.sbr_app.oanalysis.InsertSetup(
            "HfssDriven",
            [
                f"NAME:{SETUP_NAME}",
                "IsEnabled:=", True,
                ["NAME:MeshLink", "ImportMesh:=", False],
                "IsSbrRangeDoppler:=", False,
                "RayDensityPerWavelength:=", ray_density,
                "MaxNumberOfBounces:=", max_bounces,
                "RadiationSetup:=", "",
                "PTDUTDSimulationSettings:=", "PTD Correction + UTD Rays",
                "EnableCWRays:=", True,
                ["NAME:Sweeps", sweep_args],
                "ComputeFarFields:=", False,
            ],
        )  # fmt: skip

    def check_target_far_field(self, target_type):
        """Warns when the target range is shorter than the target's own far-field distance."""
        if target_type.lower() == "none":
            return
        sweep = self._source_sweep()
        start = _parse_frequency_ghz(sweep["RangeStart"])
        stop = _parse_frequency_ghz(sweep["RangeEnd"])
        centre_ghz = (start + stop) / 2.0 if start and stop else DEFAULT_CENTRE_FREQUENCY_GHZ
        wavelength_mm = SPEED_OF_LIGHT_MM_S / (centre_ghz * 1e9)
        size_mm = _length_mm(self.sbr_app, self.sbr_app["targetSize"])
        # AEDT spells metres 'meter'; _length_mm parses 'm'
        distance = re.sub(r"meter$", "m", self.sbr_app["targetDistance"].strip())
        distance_mm = _length_mm(self.sbr_app, distance)
        aperture_mm = {
            "sphere": 2.0 * size_mm,
            "dihedral": math.sqrt(3.0) * size_mm,
            "trihedral": math.sqrt(2.0) * size_mm,
        }[target_type.lower()]
        far_field_mm = 2.0 * aperture_mm**2 / wavelength_mm
        if distance_mm < far_field_mm:
            print(
                f"  Warning: target range {distance_mm / 1e3:.2f} m is inside the {target_type}'s"
                f" far-field distance 2D^2/lambda = {far_field_mm / 1e3:.2f} m"
                f" at {centre_ghz:.1f} GHz."
            )

    def create_reports(self):
        """
        (g) Quick-look reports: |S| in dB for every Tx/Rx pair (total and scattered), and per
        dual-polarised Tx/Rx element pair the co-pol ratio S_VV/S_HH (dB and phase) and the
        cross-pol isolation of the scattered (target-only) solution.
        """
        total = f"{SETUP_NAME} : {SWEEP_NAME}"
        scattered = f"{SETUP_NAME} : {SCATTERED_SWEEP_NAME}"
        pairs = [f"S({rx},{tx})" for tx in self.tx_ports for rx in self.rx_ports]
        reports = {
            REPORT_S_DB: (total, [f"dB({s})" for s in pairs]),
            REPORT_S_DB_SCATTERED: (scattered, [f"dB({s})" for s in pairs]),
        }

        imbalance, isolation = [], []
        for tx in self._dual_pol_elements(self.tx_ports).values():
            for rx in self._dual_pol_elements(self.rx_ports).values():
                s_hh, s_vv = f"S({rx['H']},{tx['H']})", f"S({rx['V']},{tx['V']})"
                s_hv, s_vh = f"S({rx['V']},{tx['H']})", f"S({rx['H']},{tx['V']})"
                imbalance += [f"dB({s_vv}/{s_hh})", f"ang_deg({s_vv}/{s_hh})"]
                isolation += [f"dB({s_hv}/{s_hh})", f"dB({s_vh}/{s_vv})"]
        if imbalance:
            reports[REPORT_COPOL_IMBALANCE] = (scattered, imbalance)
            reports[REPORT_CROSSPOL_ISOLATION] = (scattered, isolation)

        # Native CreateReport: PyAEDT's create_report rejects '_Scattered' before the first solve
        report_setup = self.sbr_app.odesign.GetModule("ReportSetup")
        families = ["Freq:=", ["All"]]
        for variable in self.sbr_app.variable_manager.independent_design_variable_names:
            families += [f"{variable}:=", ["Nominal"]]
        for name, (solution, expressions) in reports.items():
            print(f"Creating report '{name}' ({solution}, {len(expressions)} traces)...")
            try:
                report_setup.CreateReport(
                    name,
                    "Modal Solution Data",
                    "Rectangular Plot",
                    solution,
                    ["Domain:=", "Sweep"],
                    families,
                    ["X Component:=", "Freq", "Y Component:=", expressions],
                )
            except Exception as e:
                print(f"  Warning: failed to create report '{name}' ({e}).")

    def run_solve(self, simulate=None):
        """(h) Solves the SBR+ setup (and the linked HFSS design where its solution is missing)."""
        question = (
            f"Start the SBR+ simulation of '{self.sbr_design_name}' now? This also solves "
            f"'{self.source_design_name}' if its linked solution is missing"
        )
        if not _confirm(question, preset=simulate, default=False):
            print("Simulation not started.")
            return
        print(f"Solving SBR+ design '{self.sbr_design_name}'...")
        self.sbr_app.analyze_setup(SETUP_NAME)
        self.sbr_app.save_project()

    # ------------------------------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------------------------------

    def _default_source_solution(self):
        """First '<setup> : <sweep>' of the HFSS design, falling back to its first solution."""
        solutions = [s for s in self.source_app.existing_analysis_sweeps if " : " in s]
        sweeps = [s for s in solutions if "Adaptive" not in s and " - " not in s]
        if not sweeps and not solutions:
            raise ValueError(f"'{self.source_design_name}' has no analysis setup to link.")
        return (sweeps or solutions)[0]

    def _source_sweep(self):
        """Frequency range of the linked HFSS sweep as SBR+ sweep properties."""
        setup_name, sweep_name = [s.strip() for s in self.source_solution.split(":", 1)]
        setup = self.source_app.get_setup(setup_name)
        if sweep_name not in setup.get_sweep_names():
            raise ValueError(
                f"Linked solution '{self.source_solution}' is not a frequency sweep; link a sweep "
                "so that the SBR+ frequencies can match it."
            )
        props = setup.get_sweep(sweep_name).props
        if props.get("RangeType") in ("LinearCount", "LinearStep", "LogScale"):
            ranges = props
        else:
            # Interpolating/discrete sweeps with sub-ranges keep them under 'SweepRanges'
            subranges = props.get("SweepRanges", {}).get("Subrange", [])
            ranges = subranges[0] if isinstance(subranges, list) and subranges else subranges
        if not ranges or "RangeStart" not in ranges:
            raise ValueError(f"Cannot read the frequency range of '{self.source_solution}'.")
        sweep = {key: ranges[key] for key in ("RangeType", "RangeStart", "RangeEnd")}
        for key in ("RangeCount", "RangeStep", "RangeSamples"):
            if key in ranges:
                sweep[key] = ranges[key]
        return sweep

    def _visualization_objects(self):
        """All model solids and sheets of the HFSS design except the radiation boundary objects."""
        src = self.source_app
        excluded = set(RADIATION_OBJECT_NAMES)
        for boundary in src.boundaries:
            if boundary.type in RADIATION_BOUNDARY_TYPES:
                for obj in boundary.props.get("Objects", []):
                    excluded.add(src.modeler.objects[obj].name if isinstance(obj, int) else obj)
        objects = list(src.modeler.solid_names) + list(src.modeler.sheet_names)
        return [name for name in objects if name not in excluded]

    def _dual_pol_elements(self, excitations):
        """Groups excitations into {element key: {'H': name, 'V': name}} for dual-pol elements."""
        prefix = f"{self.antenna_instance_name}_"
        elements = {}
        for excitation in excitations:
            _, pol, key = _parse_port(excitation[len(prefix) :])
            elements.setdefault(key, {})[pol] = excitation
        return {key: pols for key, pols in elements.items() if set(pols) == {"H", "V"}}

    # ------------------------------------------------------------------------------------------
    # Results
    # ------------------------------------------------------------------------------------------

    def export_results(self, output_path, sweep_variable=None, solution=None):
        """
        Writes the solved S-parameters to an .npz for MimoPolarimetryProcessor.from_npz, by
        default from the scattered (target-only) solution and for all solved values of
        `sweep_variable` (see sbr_results.export_sbr_parametric_sweep).
        """
        from pomaar.simulator.sbr_results import export_sbr_parametric_sweep

        if not self.sbr_app:
            raise RuntimeError("SBR+ design is not set up.")
        solution = solution or f"{SETUP_NAME} : {SCATTERED_SWEEP_NAME}"
        output_path = os.path.abspath(output_path)
        print(f"Exporting '{solution}' to {output_path}...")
        return export_sbr_parametric_sweep(self.sbr_app, output_path, solution, sweep_variable)

    def close(self):
        """Releases AEDT desktop session connection."""
        if self.desktop_session:
            if self.is_new_desktop:
                self.desktop_session.close_desktop()
                print("AEDT desktop session closed.")
            else:
                self.desktop_session.release_desktop(close_projects=False, close_on_exit=False)
                print("AEDT desktop session connection released.")
            self.desktop_session = None


def _with_unit(value, unit):
    """Appends a unit to bare numbers; strings with units or expressions are kept as they are."""
    if isinstance(value, (int, float)):
        return f"{value}{unit}"
    text = str(value).strip()
    return f"{text}{unit}" if re.fullmatch(r"[-+]?[\d.]+(?:[eE][-+]?\d+)?", text) else text


def _euler_zyz_cs_args(origin, angles):
    """CreateRelativeCS parameters for an Euler ZYZ CS; angles are (phi, theta, psi)."""
    phi, theta, psi = angles
    return [
        "NAME:RelativeCSParameters",
        "Mode:=", "Euler Angle ZYZ",
        "OriginX:=", origin[0],
        "OriginY:=", origin[1],
        "OriginZ:=", origin[2],
        "Psi:=", psi,
        "Theta:=", theta,
        "Phi:=", phi,
    ]  # fmt: skip


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Builds an SBR+ design linking an HFSS antenna design to a canonical target."
    )
    parser.add_argument("project_path", help="Path to the AEDT project (.aedt)")
    parser.add_argument("source_design_name", help="Solved HFSS antenna design in that project")
    parser.add_argument(
        "--design", default=None, help="SBR+ design name (default: Linked<source_design_name>)"
    )
    parser.add_argument(
        "--target", choices=TARGET_TYPES, default="dihedral", help="Target type (default: dihedral)"
    )
    parser.add_argument(
        "--size",
        default=DEFAULT_TARGET_SIZE,
        help=f"Sphere radius or dihedral/trihedral edge length (default: {DEFAULT_TARGET_SIZE})",
    )
    parser.add_argument(
        "--roughness",
        type=float,
        default=0.0,
        help="SBR+ surface roughness of the target's Perfect E boundary (default: 0)",
    )
    parser.add_argument(
        "--height-deviation",
        type=float,
        default=0.0,
        help="SBR+ surface height standard deviation in mm of that boundary (default: 0)",
    )
    parser.add_argument(
        "--distance",
        default=DEFAULT_TARGET_DISTANCE,
        help=f"Target distance from AntennaCS (default: {DEFAULT_TARGET_DISTANCE})",
    )
    parser.add_argument("--azimuth", default="0deg", help="Target azimuth (default: 0deg)")
    parser.add_argument(
        "--elevation", default="0deg", help="Target polar angle off boresight (default: 0deg)"
    )
    parser.add_argument(
        "--roll", default="0deg", help="Target roll about the line of sight (default: 0deg)"
    )
    parser.add_argument(
        "--field-type",
        choices=FIELD_TYPES,
        default=None,
        help="Linked antenna field type (default: ask; far field when not interactive)",
    )
    parser.add_argument(
        "--map-variables",
        nargs="+",
        default=[],
        metavar="VARIABLE",
        help="HFSS design variables mapped to SBR+ design variables of the same name",
    )
    parser.add_argument(
        "--source-solution",
        default=None,
        help="Linked HFSS solution, e.g. 'Setup : Sweep' (default: first frequency sweep)",
    )
    parser.add_argument("--ray-density", type=float, default=4, help="Rays per wavelength")
    parser.add_argument("--max-bounces", type=int, default=5, help="Maximum number of bounces")
    parser.add_argument(
        "--overwrite", action="store_true", default=None, help="Replace an existing SBR+ design"
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
    parser.add_argument(
        "--port",
        type=int,
        default=None,
        help="gRPC port of the AEDT session (default: the running AEDT GUI's port, else 50051)",
    )
    args = parser.parse_args(argv)

    manager = SbrSimulationManager(
        project_path=args.project_path,
        source_design_name=args.source_design_name,
        sbr_design_name=args.design,
        grpc_port=args.port,
    )
    try:
        manager.connect_desktop()
        manager.setup_sbr_design(
            target_type=args.target,
            target_size=args.size,
            target_roughness=args.roughness,
            target_height_deviation=args.height_deviation,
            target_distance=args.distance,
            target_azimuth=args.azimuth,
            target_elevation=args.elevation,
            target_roll=args.roll,
            field_type=args.field_type,
            mapped_variables=args.map_variables,
            source_solution=args.source_solution,
            ray_density=args.ray_density,
            max_bounces=args.max_bounces,
            overwrite=args.overwrite,
            simulate=args.simulate,
        )
    except BuildAbortedError as e:
        print(f"[INFO] {e}")
        return 0
    except NativeAedtRefusedError as e:
        print(f"[ERROR] {e}")
        return 2
    finally:
        manager.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
