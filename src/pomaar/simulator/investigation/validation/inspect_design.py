#!/usr/bin/env python3
"""
Read-only inspection of designs of the validation project: variables, coordinate systems,
objects and boundaries, excitations, setups and sweeps (as saved), reports and solutions.

    ansys-vnc exec python -m pomaar.simulator.investigation.validation.inspect_design DESIGN ...
"""

import argparse

from pomaar.simulator.frequency_sweeps import setup_sweeps
from pomaar.simulator.hfss_array_builder import saved_setups
from pomaar.simulator.investigation.validation.common import open_validation, release
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
    parser.add_argument("designs", nargs="+")
    args = parser.parse_args(argv)
    desktop, app = open_validation(args.designs[0])
    try:
        for design in args.designs:
            app.set_active_design(design)
            inspect(app)
    finally:
        release(desktop)


if __name__ == "__main__":
    main()
