#!/usr/bin/env python3
"""
Read-only inspection of designs in the live GUI session: variables, coordinate systems,
objects and boundaries, excitations, setups and sweeps (as saved), reports and solutions.
PROJECT is a .aedt file (opened in the GUI if needed) or the name of an open project.

    ansys-vnc exec python -m pomaar.agent.inspect_design PROJECT DESIGN ...
"""

import argparse
from pathlib import Path
import sys

from ansys.aedt.core import Hfss

from pomaar.simulator.frequency_sweeps import setup_sweeps
from pomaar.simulator.hfss_array_builder import (
    MissingDesignError,
    connect_live_session,
    require_designs,
    saved_setups,
)
from pomaar.simulator.sbr_results import scattered_solution


def inspect(app):
    print(f"\n===== {app.design_name} ({app.solution_type})")
    print("variables:", {v: app[v] for v in app.variable_manager.design_variable_names})
    for cs in app.modeler.coordinate_systems:
        props = {k: v for k, v in dict(cs.props).items() if k not in ("Attributes",)}
        print(f"CS {cs.name} (ref {cs.ref_cs}): {props}")
    print("objects:", sorted(app.modeler.object_names))
    print("boundaries:", sorted(b.name for b in app.boundaries))
    print("excitations:", sorted(app.excitation_names))
    for setup, props in saved_setups(app).items():
        sweeps = {
            name: {k: v for k, v in s.items() if "Range" in k or k in ("Type", "SaveRadFields")}
            for name, s in setup_sweeps(props).items()
        }
        print(f"setup {setup}: sweeps {sweeps}")
    print("reports:", sorted(app.post.all_report_names))
    print("solutions:", list(app.existing_analysis_sweeps))
    if "SBR" in str(app.solution_type):
        try:
            print("scattered solution:", scattered_solution(app))
        except Exception as error:  # noqa: BLE001 (a report, not a check)
            print("scattered solution: none", error)
        try:
            tx_rx = app.oboundary.GetSBRTxRxSettings()
            print("SBR Tx/Rx:", tx_rx)
        except Exception as error:  # noqa: BLE001
            print("SBR Tx/Rx: unreadable", error)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("project", help=".aedt file or the name of an open project")
    parser.add_argument("designs", nargs="+")
    args = parser.parse_args(argv)
    desktop = connect_live_session()
    try:
        # Checked first: opening a missing design by name would create it
        try:
            require_designs(desktop, args.project, args.designs)
        except MissingDesignError as error:
            sys.exit(f"[ERROR] {error}")
        name = Path(args.project).stem if args.project.endswith(".aedt") else args.project
        app = Hfss(project=name, design=args.designs[0], new_desktop=False, close_on_exit=False)
        for design in args.designs:
            app.set_active_design(design)
            inspect(app)
    finally:
        desktop.release_desktop(close_projects=False, close_on_exit=False)


if __name__ == "__main__":
    main()
