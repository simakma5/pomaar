#!/usr/bin/env python3
"""
SBR+ result export: pulls the swept S-parameters of a solved (parametric) SBR+ design into an
.npz file, so that all further processing runs offline in `polarimetry_processor`.

Run inside the AEDT container, attached to the running GUI session:

    ansys-vnc exec sbr_results <project> <design> <out.npz> \
        [--setup "Setup : Sweep_Scattered"] [--sweep-variable copolarSpacingLambda] \
        [--port-table <ArrayDesign>_ports.json]

The scattered (target-only) solution is the default; see the README on why it is the
simulation's counterpart of range gating out the antenna coupling in a measurement.
"""

import argparse
import glob
import json
import os
import sys

from ansys.aedt.core import Hfss
import numpy as np

from pomaar.simulator.array_templates import attach_port_table, s_parameter_ports
from pomaar.simulator.hfss_array_builder import (
    NativeAedtRefusedError,
    NoLiveSessionError,
    connect_live_session,
)

SCATTERED_SOLUTION = "Setup : Sweep_Scattered"


def find_port_table(expressions, search_dir, ports_path=None):
    """
    The builder port table ('<ArrayDesign>_ports.json', written next to the project) that matches
    the excitations of `expressions`, as (rows with their 'excitation', path), or (None, None).
    `ports_path` names the table; otherwise every '*_ports.json' in `search_dir` is tried.
    """
    excitations = [p for e in expressions if (pair := s_parameter_ports(e)) for p in pair]
    if ports_path:
        candidates = [ports_path]
    else:
        candidates = sorted(glob.glob(os.path.join(search_dir, "*_ports.json")))
    for path in candidates:
        try:
            with open(path, encoding="utf-8") as f:
                rows = attach_port_table(json.load(f), excitations)
        except (OSError, ValueError, KeyError) as e:
            print(f"  Warning: cannot read port table {path} ({e})")
            continue
        if rows:
            return rows, path
    if ports_path:
        raise ValueError(f"Port table {ports_path} does not match the exported excitations.")
    return None, None


def export_sbr_parametric_sweep(
    app, output_path, setup_sweep=SCATTERED_SOLUTION, sweep_variable=None, ports_path=None
):
    """
    Writes every S-parameter of `setup_sweep` for all solved values of `sweep_variable`.

    A solved SBR+ design offers 'Setup : Sweep' (total fields), 'Setup : Sweep_Incident' (direct
    Tx->Rx coupling only) and 'Setup : Sweep_Scattered' (target only, the default); target
    signatures need the scattered one.

    S is read as AEDT's real and imaginary parts and stored complex: that is AEDT's own form,
    exact, and free of the -inf dB of a zero channel and of phase wrapping (magnitude/phase data
    from elsewhere rebuild with `polarimetry_processor.complex_from_db_deg`).

    The .npz holds `s` (n_variations, n_expressions, n_frequencies) complex, `frequencies_hz`,
    `sweep_values`, `expressions` (AEDT names, e.g. 'S(HornCluster1_RxH,HornCluster1_TxV)', i.e.
    S(receive, transmit)), and `metadata` (JSON: project, design, setup, sweep variable and the
    nominal design variables). Without a sweep variable only the nominal variation is written.
    When the linked antenna was built by `hfss_array_builder`, its port table (`ports_path`, or
    found next to the project) is added to the metadata as 'ports': per excitation the builder
    port, role, polarization, field sign and phase centre (mm, at build time), which
    `MimoPolarimetryProcessor` uses to correct the signs and as virtual-channel positions.
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
    rows, table_path = find_port_table(expressions, app.project_path, ports_path)
    if rows:
        metadata["ports"] = rows
        metadata["port_table"] = table_path
        print(f"Attached the builder port table {table_path} ({len(rows)} ports)")
    else:
        print("No builder port table found; port signs and positions are not attached")
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
    parser.add_argument(
        "--setup",
        default=SCATTERED_SOLUTION,
        help=f"Setup and sweep name (default: '{SCATTERED_SOLUTION}', the target only)",
    )
    parser.add_argument("--sweep-variable", default=None, help="Parametric variable to export")
    parser.add_argument(
        "--port-table",
        default=None,
        help=(
            "Builder port table of the linked array (default: the matching '*_ports.json' next "
            "to the project)"
        ),
    )
    args = parser.parse_args()

    try:
        desktop = connect_live_session()
    except (NativeAedtRefusedError, NoLiveSessionError) as e:
        sys.exit(f"[ERROR] {e}")
    try:
        app = Hfss(project=args.project, design=args.design, new_desktop=False, close_on_exit=False)
        export_sbr_parametric_sweep(
            app, args.output, args.setup, args.sweep_variable, ports_path=args.port_table
        )
    finally:
        desktop.release_desktop(close_projects=False, close_on_exit=False)


if __name__ == "__main__":
    main()
