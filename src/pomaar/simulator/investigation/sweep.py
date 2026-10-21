#!/usr/bin/env python3
"""
Creates (or replaces) an Optimetrics parametric sweep of an investigation design and optionally
solves it, e.g.

    ansys-vnc exec python src/pomaar/simulator/investigation/sweep.py MultiRxDihedral \
        --sweep bistaticAngle "LIN 0.25deg 8deg 0.25deg" --sweep rxAim "LIN 0 1 1" \
        --sweep targetRoll "LIN 0deg 90deg 22.5deg" --solve

Every --sweep is one AEDT sweep definition (Data: 'LIN start stop step', 'LINC start stop
count' or a single value); the variations are their Cartesian product.
"""

import argparse
import time

from pomaar.simulator.investigation.project import add_project_argument, open_project, release
from pomaar.simulator.sbr_results import scattered_solution

SWEEP_NAME = "BistaticSweep"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("design")
    parser.add_argument("--name", default=SWEEP_NAME)
    parser.add_argument(
        "--sweep", nargs=2, action="append", required=True, metavar=("VARIABLE", "DATA")
    )
    parser.add_argument("--solve", action="store_true")
    parser.add_argument(
        "--validation", action="store_true", help="Use the pipeline-validation project instead"
    )
    add_project_argument(parser)
    args = parser.parse_args(argv)

    if args.validation:
        from pomaar.simulator.investigation.validation.common import open_validation

        desktop, app = open_validation(args.design)
    else:
        desktop, app = open_project(args.design, args.project)
    try:
        setup = scattered_solution(app).split(" : ")[0]
        optimetrics = app.odesign.GetModule("Optimetrics")
        if args.name in list(optimetrics.GetSetupNames()):
            optimetrics.DeleteSetups([args.name])
            print(f"Replaced the parametric setup '{args.name}'.")
        definitions = [
            [
                "NAME:SweepDefinition",
                "Variable:=", variable,
                "Data:=", data,
                "OffsetF1:=", False,
                "Synchronize:=", 0,
            ]
            for variable, data in args.sweep
        ]  # fmt: skip
        optimetrics.InsertSetup(
            "OptiParametric",
            [
                f"NAME:{args.name}",
                "IsEnabled:=", True,
                [
                    "NAME:ProdOptiSetupDataV2",
                    "SaveFields:=", False,
                    "CopyMesh:=", False,
                    "SolveWithCopiedMeshOnly:=", True,
                ],
                ["NAME:StartingPoint"],
                "Sim. Setups:=", [setup],
                ["NAME:Sweeps"] + definitions,
                ["NAME:Sweep Operations"],
                ["NAME:Goals"],
            ],
        )  # fmt: skip
        print(f"'{args.name}' on '{setup}': {args.sweep}")
        app.save_project()
        if args.solve:
            start = time.time()
            app.analyze_setup(args.name)
            print(f"Solved '{args.name}' in {(time.time() - start) / 60:.1f} min")
            app.save_project()
    finally:
        release(desktop)


if __name__ == "__main__":
    main()
