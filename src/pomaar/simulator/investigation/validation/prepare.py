#!/usr/bin/env python3
"""V0: cuts the unit cell's sweep in the validation copy to a few points (fast array solves)."""

import argparse

from pomaar.simulator.hfss_array_builder import saved_setups
from pomaar.simulator.investigation.validation.common import open_validation, release

parser = argparse.ArgumentParser()
parser.add_argument("--design", default="ApertureCoupledPatch")
parser.add_argument("--points", type=int, default=5)
args = parser.parse_args()

desktop, app = open_validation(args.design)
try:
    setups = saved_setups(app)
    setup = next(iter(setups))
    sweep = list(app.oanalysis.GetSweeps(setup))[0]
    props = setups[setup]["Sweeps"][sweep] if "Sweeps" in setups[setup] else None
    print("setup", setup, "sweep", sweep, "cs", [c.name for c in app.modeler.coordinate_systems])
    sweep_type = (props or {}).get("Type", "Interpolating")
    app.oanalysis.EditFrequencySweep(
        setup,
        sweep,
        [
            f"NAME:{sweep}",
            "IsEnabled:=", True,
            "RangeType:=", "LinearCount",
            "RangeStart:=", "77GHz",
            "RangeEnd:=", "81GHz",
            "RangeCount:=", args.points,
            "Type:=", sweep_type,
            "SaveFields:=", False,
            "SaveRadFields:=", False,
        ],
    )  # fmt: skip
    after = saved_setups(app)[setup]
    edited = after["Sweeps"][sweep]
    print("sweep after edit:", {k: v for k, v in edited.items() if "Range" in k or k == "Type"})
    app.save_project()
finally:
    release(desktop)
