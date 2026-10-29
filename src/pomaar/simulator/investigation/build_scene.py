#!/usr/bin/env python3
"""
Builds a multi-receiver bistaticity scene as a new SBR+ design of the investigation project
(x boresight, ISO 8855: x forward, y left, z up; H along y, V along z):

* one dual-polarized Tx site at the origin and `--rx-count` dual-polarized Rx sites, all
  parametric beams (`verticalBeamWidth` x `horizontalBeamWidth`);
* Rx site k at y_k = targetRange tan(k bistaticAngle), so that a boresight target sees it at
  the bistatic angle k bistaticAngle (k = 0: monostatic, co-located with the Tx site);
* each Rx turned by yaw -rxAim k bistaticAngle (rxAim = 1: aimed at a boresight target,
  0: on boresight like a real array);
* H at z = +spacingHV/2 and V at z = -spacingHV/2 in every site;
* TargetCS and the canonical target as `sbr_simulator` builds them (or a thin PEC strip along
  H for the polarization-label check, or a square PEC plate of edge targetSize facing the
  radar), and an SBR+ setup like the original design's.

    ansys-vnc exec python src/pomaar/simulator/investigation/build_scene.py MultiRxTrihedral \
        --target trihedral [--rx-count 3] [--bistatic-angle 1deg] [--overwrite]

`--layout pairs` builds the split-pair scene for the differential reflectivity instead: four
single-polarized antennas on the y axis in the order TxV, RxV, (gap), TxH, RxH,

    TxV at y = g/2 + b,  RxV at g/2,  TxH at -g/2,  RxH at -g/2 - b,

with the pair base b = 2 targetRange tan(pairAngle / 2) and the gap g = gapRatio b, so that VV
and HH are each measured by one pair of base b (bistatic angle about pairAngle) whose bisectors
lie (b + g) / targetRange apart. The VV pair is the mirror image of the HH pair about the x-z
plane, so both have the same path length and the same beam gains for a boresight target. A
co-located dual-polarized reference pair of the same base (Tx_R at y = b/2, Rx_R at -b/2)
measures the full matrix with its bisector on boresight. All antennas look along boresight.
"""

import argparse
import sys

from ansys.aedt.core import Hfss
from ansys.aedt.core.generic.constants import Plane

from pomaar.simulator.frequency_sweeps import native_sweep_arguments
from pomaar.simulator.investigation.project import (
    add_project_argument,
    open_project,
    release,
    reopen_project,
)
from pomaar.simulator.sbr_simulator import (
    TARGET_CS,
    TARGET_TYPES,
    SbrSimulationManager,
    _euler_zyz_cs_args,
)

TARGETS = tuple(t for t in TARGET_TYPES if t != "none") + ("strip", "plate", "none")
POLARIZATIONS = {"H": "Horizontal", "V": "Vertical"}
REFERENCE_DESIGN = "BistaticityAnalysis"


def site_names(rx_count):
    return ["Tx"] + [f"Rx_{k}" for k in range(rx_count)]


def replace_design(desktop, project_app, name, overwrite):
    """
    (Re)creates the SBR+ design `name` with native calls (see the CLAUDE.md pitfalls). A
    deleted design's result files stay open in AEDT, so the project is reopened before the
    name is used again (`reopen_project`).
    """
    oproject = project_app.oproject
    if name in oproject.GetTopDesignList():
        if not overwrite:
            sys.exit(f"[ERROR] Design '{name}' exists; pass --overwrite to rebuild it.")
        retired = f"{name}Deleted"
        oproject.GetDesign(name).RenameDesignInstance(name, retired)
        oproject.DeleteDesign(retired)
        print(f"Deleted the existing design '{name}'; reopening the project...")
        project_app = reopen_project(desktop, project_app, REFERENCE_DESIGN)
        oproject = project_app.oproject
    oproject.InsertDesign("HFSS", name, "SBR+", "")
    if name not in oproject.GetTopDesignList():
        raise RuntimeError(f"AEDT did not create design '{name}'.")
    app = Hfss(
        project=project_app.project_name, design=name, new_desktop=False, close_on_exit=False
    )
    app.modeler.model_units = "mm"
    return app


def create_strip(app, length="targetSize", width="0.5mm", name="TargetStrip"):
    """
    A PEC rectangle in TargetCS, centred, along its local x axis (H at zero roll) and facing
    the radar: a thin strip by default, a square plate with width = length.
    """
    modeler = app.modeler
    modeler.set_working_coordinate_system(TARGET_CS)
    strip = modeler.create_rectangle(
        Plane.XY, [f"-{length}/2", f"-({width})/2", "0mm"], [length, width], name=name
    )
    modeler.set_working_coordinate_system("Global")
    app.assign_perfect_e(strip.name, name=f"{name}PerfectE")


def create_antennas(app, rx_count):
    """Antenna CSs and parametric beams of every site; returns ({site: {pol: instance}})."""
    editor = app.modeler.oeditor
    sites = {}
    for k, site in enumerate(site_names(rx_count)):
        index = k - 1  # Rx index; -1 for the Tx site
        y = "0mm" if index <= 0 else f"targetRange*tan({index}*bistaticAngle)"
        yaw = "0deg" if index <= 0 else f"-rxAim*{index}*bistaticAngle"
        sites[site] = {}
        for pol, z in (("H", "spacingHV/2"), ("V", "-spacingHV/2")):
            cs = f"{site}_{pol}_CS"
            app.modeler.set_working_coordinate_system("Global")
            editor.CreateRelativeCS(
                _euler_zyz_cs_args(["0mm", y, z], [yaw, "0deg", "0deg"]),
                ["NAME:Attributes", "Name:=", cs],
            )
            app.modeler.set_working_coordinate_system("Global")
            antenna = app.create_sbr_antenna(
                app.SbrAntennas.ParametricBeam,
                target_cs=cs,
                parameters={
                    "Size": "0.1mm",
                    "Polarization": POLARIZATIONS[pol],
                    "Vertical BeamWidth": "verticalBeamWidth",
                    "Horizontal BeamWidth": "horizontalBeamWidth",
                },
                name=f"{site}_{pol}",
            )
            sites[site][pol] = antenna.name
    return sites


def create_pair_antennas(app):
    """
    The split pairs (Tx_A_V, Rx_A_V; Tx_B_H, Rx_B_H) and the co-located reference pair (Tx_R
    and Rx_R, H and V each) of `--layout pairs`, all on boresight; returns {instance: y}.
    """
    editor = app.modeler.oeditor
    base = "(2*targetRange*tan(pairAngle/2))"
    gap = f"(gapRatio*{base})"
    placements = {
        "Tx_A_V": f"{gap}/2+{base}",
        "Rx_A_V": f"{gap}/2",
        "Tx_B_H": f"-{gap}/2",
        "Rx_B_H": f"-{gap}/2-{base}",
        "Tx_R_H": f"{base}/2",
        "Tx_R_V": f"{base}/2",
        "Rx_R_H": f"-{base}/2",
        "Rx_R_V": f"-{base}/2",
    }
    for name, y in placements.items():
        cs = f"{name}_CS"
        app.modeler.set_working_coordinate_system("Global")
        editor.CreateRelativeCS(
            _euler_zyz_cs_args(["0mm", y, "0mm"], ["0deg", "0deg", "0deg"]),
            ["NAME:Attributes", "Name:=", cs],
        )
        app.modeler.set_working_coordinate_system("Global")
        app.create_sbr_antenna(
            app.SbrAntennas.ParametricBeam,
            target_cs=cs,
            parameters={
                "Size": "0.1mm",
                "Polarization": POLARIZATIONS[name[-1]],
                "Vertical BeamWidth": "verticalBeamWidth",
                "Horizontal BeamWidth": "horizontalBeamWidth",
            },
            name=name,
        )
    return placements


def assign_tx_rx(app):
    """Every Tx port transmits to every Rx port; returns (tx ports, rx ports)."""
    # SetSBRTxRxSettings silently ignores PyAEDT's mode suffix (':1')
    ports = [name.split(":")[0] for name in app.excitation_names]
    tx = [p for p in ports if p.startswith("Tx_")]
    rx = [p for p in ports if p.startswith("Rx_")]
    args = ["NAME:SBRTxRxSettings"]
    for index, port in enumerate(tx):
        args.append(
            [f"NAME:Tx/Rx List {index}", "Tx Antenna:=", port, "Rx Antennas:=", ",".join(rx)]
        )
    app.oboundary.SetSBRTxRxSettings(args)
    return tx, rx


def create_setup(app, start, stop, step, ray_density, max_bounces):
    """The original design's SBR+ setup; returns '<setup> : <sweep>' as AEDT named them."""
    ranges = [{"RangeType": "LinearStep", "RangeStart": start, "RangeEnd": stop, "RangeStep": step}]
    analysis = app.oanalysis
    existing = set(analysis.GetSetups())
    analysis.InsertSetup(
        "HfssDriven",
        [
            "NAME:Setup",
            "IsEnabled:=", True,
            ["NAME:MeshLink", "ImportMesh:=", False],
            "IsSbrRangeDoppler:=", False,
            "RayDensityPerWavelength:=", ray_density,
            "MaxNumberOfBounces:=", max_bounces,
            "RadiationSetup:=", "",
            "PTDUTDSimulationSettings:=", "PTD Correction + UTD Rays",
            ["NAME:Sweeps", native_sweep_arguments("Sweep", ranges)],
            "ComputeFarFields:=", False,
        ],
    )  # fmt: skip
    setup = sorted(set(analysis.GetSetups()) - existing)[0]
    return f"{setup} : {list(analysis.GetSweeps(setup))[0]}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("design", help="Name of the SBR+ design to create")
    parser.add_argument("--target", choices=TARGETS, default="trihedral")
    parser.add_argument("--size", default="40mm", help="Target size (sbr_simulator's meaning)")
    parser.add_argument("--range", default="4m", help="targetRange (default: 4m)")
    parser.add_argument("--roll", default="0deg", help="targetRoll (default: 0deg)")
    parser.add_argument(
        "--layout",
        choices=("ladder", "pairs"),
        default="ladder",
        help="Receiver ladder (default) or the split pairs for the differential reflectivity",
    )
    parser.add_argument("--rx-count", type=int, default=3, help="Rx sites (default: 3)")
    parser.add_argument("--pair-angle", default="1deg", help="pairAngle (pairs; default: 1deg)")
    parser.add_argument("--gap-ratio", default="0", help="gapRatio (pairs; default: 0)")
    parser.add_argument("--bistatic-angle", default="1deg", help="bistaticAngle (default: 1deg)")
    parser.add_argument("--rx-aim", type=int, choices=(0, 1), default=1, help="rxAim (default: 1)")
    parser.add_argument("--spacing-hv", default="0mm", help="spacingHV (default: 0mm, dual-pol)")
    parser.add_argument("--beam-width", default="60deg", help="Both beam widths (default: 60deg)")
    parser.add_argument("--start", default="74GHz")
    parser.add_argument("--stop", default="84GHz")
    parser.add_argument("--step", default="0.01GHz")
    parser.add_argument("--ray-density", type=float, default=4)
    parser.add_argument("--max-bounces", type=int, default=5)
    parser.add_argument("--overwrite", action="store_true")
    add_project_argument(parser)
    args = parser.parse_args(argv)

    desktop, reference = open_project(REFERENCE_DESIGN, args.project)
    try:
        app = replace_design(desktop, reference, args.design, args.overwrite)
        if args.layout == "pairs":
            variables = (("pairAngle", args.pair_angle), ("gapRatio", args.gap_ratio))
        else:
            variables = (
                ("spacingHV", args.spacing_hv),
                ("bistaticAngle", args.bistatic_angle),
                ("rxAim", str(args.rx_aim)),
            )
        variables = (
            ("verticalBeamWidth", args.beam_width),
            ("horizontalBeamWidth", args.beam_width),
        ) + variables
        for name, value in variables:
            app[name] = value

        manager = SbrSimulationManager(
            str(app.project_file), REFERENCE_DESIGN, args.design, boresight="x"
        )
        manager.sbr_app = app
        manager.create_coordinate_systems(
            target_range=args.range, target_roll=args.roll, create_antenna_cs=False
        )
        if args.target in ("strip", "plate"):
            app["targetSize"] = args.size
            if args.target == "strip":
                create_strip(app)
            else:
                create_strip(app, width="targetSize", name="TargetPlate")
        else:
            manager.create_target(target_type=args.target, size=args.size)

        if args.layout == "pairs":
            sites = create_pair_antennas(app)
        else:
            sites = create_antennas(app, args.rx_count)
        tx, rx = assign_tx_rx(app)
        solution = create_setup(
            app, args.start, args.stop, args.step, args.ray_density, args.max_bounces
        )
        app.modeler.set_working_coordinate_system("Global")
        app.save_project()
        print(f"Built '{args.design}': sites {sites}")
        print(f"Tx {tx} -> Rx {rx}; setup '{solution}'")
    finally:
        release(desktop)


if __name__ == "__main__":
    main()
