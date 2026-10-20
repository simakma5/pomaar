#!/usr/bin/env python3
"""
Solves a design of the sandbox project (nominal variation, or an existing parametric sweep)
and exports its scattered solution to <sandbox>/exports/<name>.npz:

    ansys-vnc exec python src/pomaar/simulator/investigation/solve.py MultiRxTrihedral \
        --export trihedral_nominal [--set bistaticAngle=10deg ...] [--ray-density 8] \
        [--ptd-utd None] \
        [--parametric NAME --sweep-variable VAR]
"""

import argparse
import time

from pomaar.simulator.investigation.sandbox import SANDBOX_DIR, open_sandbox, release
from pomaar.simulator.sbr_results import export_sbr_parametric_sweep, scattered_solution


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("design")
    parser.add_argument("--export", required=True, help="Export name (no extension)")
    parser.add_argument("--set", nargs="*", default=[], metavar="VAR=VALUE")
    parser.add_argument("--ray-density", type=float, default=None)
    parser.add_argument(
        "--ptd-utd",
        choices=("None", "PTD Correction", "PTD Correction + UTD Rays"),
        default=None,
        help="Edge-diffraction model of the setup (PTDUTDSimulationSettings)",
    )
    parser.add_argument("--parametric", default=None, help="Optimetrics setup to solve instead")
    parser.add_argument("--sweep-variable", default=None)
    parser.add_argument("--no-solve", action="store_true", help="Only export")
    args = parser.parse_args(argv)

    desktop, app = open_sandbox(args.design)
    try:
        for assignment in args.set:
            name, value = assignment.split("=", 1)
            app[name] = value
            print(f"{name} = {app[name]}")
        setup_name = scattered_solution(app).split(" : ")[0]
        if args.ray_density is not None:
            setup = app.get_setup(setup_name)
            setup.props["RayDensityPerWavelength"] = args.ray_density
            setup.update()
            print(f"Ray density {args.ray_density} per wavelength")
        if args.ptd_utd is not None:
            setup = app.get_setup(setup_name)
            setup.props["PTDUTDSimulationSettings"] = args.ptd_utd
            setup.update()
            print(f"Edge diffraction: {args.ptd_utd}")
        if not args.no_solve:
            start = time.time()
            if args.parametric:
                app.analyze_setup(args.parametric)
            else:
                app.analyze_setup(setup_name)
            print(f"Solved in {time.time() - start:.0f} s")
        changed = args.ray_density is not None or args.ptd_utd is not None
        if args.set or changed or not args.no_solve:
            app.save_project()
        exports = SANDBOX_DIR / "exports"
        exports.mkdir(exist_ok=True)
        export_sbr_parametric_sweep(
            app, str(exports / f"{args.export}.npz"), sweep_variable=args.sweep_variable
        )
    finally:
        release(desktop)


if __name__ == "__main__":
    main()
