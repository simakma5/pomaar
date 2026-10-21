#!/usr/bin/env python3
"""B4: setup read-back. Renames ValCluster's setup, then a fresh builder finds it by itself."""

from pathlib import Path

from pomaar.simulator.hfss_array_builder import MimoHfssBuilder
from pomaar.simulator.investigation.validation.common import (
    PROJECT_PATH,
    VALIDATION_DIR,
    open_validation,
    release,
)

LAYOUT = Path(__file__).parents[2] / "layouts" / "ApertureCoupledPatch_polarimetric_cluster.yaml"

desktop, app = open_validation("ValCluster")
try:
    if "Setup" in app.oanalysis.GetSetups():
        app.oanalysis.RenameSetup("Setup", "ArraySetup")
        app.save_project()
    print("setups now:", list(app.oanalysis.GetSetups()))
    print("solutions:", [s for s in app.existing_analysis_sweeps if " - " not in s])
finally:
    release(desktop)

builder = MimoHfssBuilder(str(PROJECT_PATH), str(LAYOUT), "ValCluster")
builder.connect_desktop()
try:
    builder.target_design_app = builder._open_design("ValCluster")
    print("array_solution():", builder.array_solution())
    builder.create_post_processing_reports()
    out = VALIDATION_DIR / "b4_ValCluster.s4p"
    builder.export_coupling_s_parameters(str(out))
    print("touchstone:", out, out.is_file(), out.stat().st_size if out.is_file() else None)
finally:
    builder.close()
