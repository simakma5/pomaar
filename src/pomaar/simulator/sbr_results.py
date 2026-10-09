#!/usr/bin/env python3
"""
SBR+ result export: pulls the swept S-parameters of a solved (parametric) SBR+ design into an
.npz file, so that all further processing runs offline in `polarimetry_processor`.

Run inside the AEDT container, attached to the running GUI session:

    ansys-vnc exec .venv/bin/python -m pomaar.simulator.sbr_results <project> <design> <out.npz> \
        [--setup "Setup : Sweep_Scattered"] [--sweep-variable copolarSpacingLambda] [--port 50051]
"""

import argparse
import json
import sys

from ansys.aedt.core import Hfss
from ansys.aedt.core.generic.general_methods import grpc_active_sessions
import numpy as np


def export_sbr_parametric_sweep(app, output_path, setup_sweep="Setup : Sweep", sweep_variable=None):
    """
    Writes every S-parameter of `setup_sweep` for all solved values of `sweep_variable`.

    A solved SBR+ design offers 'Setup : Sweep' (total fields), 'Setup : Sweep_Incident' (direct
    Tx->Rx coupling only) and 'Setup : Sweep_Scattered' (target only); use the scattered one for
    target signatures.

    The .npz holds `s` (n_variations, n_expressions, n_frequencies) complex, `frequencies_hz`,
    `sweep_values`, `expressions` (AEDT names, e.g. 'S(HornCluster1_RxH,HornCluster1_TxV)', i.e.
    S(receive, transmit)), and `metadata` (JSON: project, design, setup, sweep variable and the
    nominal design variables). Without a sweep variable only the nominal variation is written.
    The result data must exist: an SBR+ design whose results were invalidated (e.g. by an edit
    of the linked source) reports "No Data Available" until it is re-solved.
    """
    expressions = app.post.available_report_quantities(
        solution=setup_sweep, quantities_category="S Parameter"
    )
    variations = {"Freq": ["All"]}
    if sweep_variable:
        variations[sweep_variable] = ["All"]
    solution_data = app.post.get_solution_data(
        expressions=expressions,
        setup_sweep_name=setup_sweep,
        variations=variations,
        primary_sweep_variable="Freq",
        report_category="Modal Solution Data",
    )
    if not solution_data:
        raise RuntimeError(f"No solution data for '{setup_sweep}'; is the design solved?")

    sweep_values = (
        list(solution_data.variation_values(sweep_variable)) if sweep_variable else [np.nan]
    )
    frequencies = np.asarray(solution_data.primary_sweep_values, dtype=float)
    s = np.zeros((len(sweep_values), len(expressions), len(frequencies)), dtype=complex)
    for variation_index, value in enumerate(sweep_values):
        solution_data.set_active_variation(variation_index)
        if sweep_variable:
            active = solution_data.active_variation[sweep_variable]
            if not np.isclose(active, value):
                raise RuntimeError(f"Variation order mismatch: {active} != {value}")
        for expression_index, expression in enumerate(expressions):
            _, real = solution_data.get_expression_data(expression, formula="real")
            _, imag = solution_data.get_expression_data(expression, formula="imag")
            s[variation_index, expression_index] = np.asarray(real) + 1j * np.asarray(imag)

    # PyAEDT returns the primary sweep in the report's unit (GHz for SBR+ sweeps)
    frequency_unit = solution_data.units_sweeps.get("Freq") or "GHz"
    scale = {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9, "THz": 1e12}[frequency_unit]
    metadata = {
        "project": app.project_name,
        "design": app.design_name,
        "setup_sweep": setup_sweep,
        "sweep_variable": sweep_variable,
        "nominal_variables": app.available_variations.nominal_values,
    }
    np.savez(
        output_path,
        s=s,
        frequencies_hz=frequencies * scale,
        sweep_values=np.asarray(sweep_values, dtype=float),
        expressions=np.asarray(expressions),
        metadata=json.dumps(metadata),
    )
    print(f"Wrote {s.shape} S-parameters ({len(expressions)} channels) to {output_path}")
    return s, frequencies * scale, sweep_values, expressions


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("project", help="Project name (open in the GUI) or .aedt path")
    parser.add_argument("design", help="Solved SBR+ design")
    parser.add_argument("output", help="Output .npz path")
    parser.add_argument("--setup", default="Setup : Sweep", help="Setup and sweep name")
    parser.add_argument("--sweep-variable", default=None, help="Parametric variable to export")
    parser.add_argument("--port", type=int, default=50051, help="gRPC port of the GUI session")
    args = parser.parse_args()

    # Only ever attach: PyAEDT would otherwise start a fresh AEDT on a free port
    if args.port not in grpc_active_sessions():
        sys.exit(f"No AEDT listening on gRPC port {args.port}; refusing to start one.")
    app = Hfss(
        project=args.project,
        design=args.design,
        new_desktop=False,
        port=args.port,
        close_on_exit=False,
    )
    try:
        export_sbr_parametric_sweep(app, args.output, args.setup, args.sweep_variable)
    finally:
        app.release_desktop(close_projects=False, close_desktop=False)


if __name__ == "__main__":
    main()
