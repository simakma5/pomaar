#!/usr/bin/env python3
"""
SBR+ result export: pulls the swept S-parameters of a solved (parametric) SBR+ design into an
.npz file, so that all further processing runs offline in `polarimetry_processor`.

Run inside the AEDT container, attached to the running GUI session:

    ansys-vnc exec sbr_results <project> <design> <out.npz> \
        [--setup "<setup> : <sweep>_Scattered"] [--sweep-variable copolarSpacingLambda] \
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
    MissingDesignError,
    NativeAedtRefusedError,
    NoLiveSessionError,
    connect_live_session,
    require_designs,
)
from pomaar.simulator.polarimetry_processor import (
    IncompleteSweepError,
    check_complete_sweep,
    identical_variations,
)
from pomaar.simulator.sbr_simulator import _degrees, antenna_cs_boresight
from pomaar.simulator.scene_frames import target_steering_angles

SCATTERED_SUFFIX = "_Scattered"


def scattered_solution(app):
    """
    The design's scattered (target-only) solution '<setup> : <sweep>_Scattered', with the names
    AEDT gave the setup and sweep (they may be de-duplicated, e.g. 'Setup1'): the design's only
    setup and its only sweep, or the first ones with a warning.
    """
    analysis = app.oanalysis
    setups = list(analysis.GetSetups())
    if not setups:
        raise RuntimeError(f"Design '{app.design_name}' has no analysis setup.")
    sweeps = list(analysis.GetSweeps(setups[0]))
    if not sweeps:
        raise RuntimeError(f"Setup '{setups[0]}' of '{app.design_name}' has no frequency sweep.")
    if len(setups) > 1 or len(sweeps) > 1:
        print(f"Warning: setups {setups}, sweeps {sweeps}; using '{setups[0]} : {sweeps[0]}'.")
    return f"{setups[0]} : {sweeps[0]}{SCATTERED_SUFFIX}"


def resolve_output_path(output_path):
    """
    The absolute .npz path of `output_path` (relative paths resolve against the current
    directory), checked before any AEDT work: its folder must exist and be writable.
    """
    path = os.path.abspath(os.path.expanduser(str(output_path)))
    if not path.endswith(".npz"):
        path += ".npz"  # np.savez would append it anyway
    folder = os.path.dirname(path)
    if not os.path.isdir(folder):
        raise FileNotFoundError(
            f"Output folder {folder} does not exist (the output path '{output_path}' is resolved "
            f"against the current directory {os.getcwd()}). Create it, or give another path."
        )
    if not os.access(folder, os.W_OK):
        raise PermissionError(f"Output folder {folder} is not writable.")
    return path


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


def scene_target(app, sweep_variable=None, sweep_values=None):
    """
    {'boresight', 'target_steering_deg'} of an `sbr_simulator` scene: the boresight read from
    AntennaCS (`sbr_simulator.antenna_cs_boresight`) and the processor's (azimuth, elevation) of
    the target per variation (`scene_frames.target_steering_angles`), so that an analysis can
    steer at the target in either frame. Empty for scenes without AntennaCS or target angles.
    """
    boresight = antenna_cs_boresight(app)
    nominal = app.available_variations.nominal_values
    if boresight is None or not {"targetAzimuth", "targetElevation"} <= set(nominal):
        return {}
    angles = []
    for value in sweep_values if sweep_variable else [None]:
        variables = dict(nominal)
        if sweep_variable in ("targetAzimuth", "targetElevation"):
            variables[sweep_variable] = value
        azimuth, elevation = (
            _degrees(variables["targetAzimuth"]),
            _degrees(variables["targetElevation"]),
        )
        if azimuth is None or elevation is None:
            return {"boresight": boresight}
        angles.append(list(target_steering_angles(azimuth, elevation, boresight)))
    return {"boresight": boresight, "target_steering_deg": angles}


def export_sbr_parametric_sweep(
    app, output_path, setup_sweep=None, sweep_variable=None, ports_path=None
):
    """
    Writes every S-parameter of `setup_sweep` for all solved values of `sweep_variable`.

    A solved SBR+ design offers '<setup> : <sweep>' (total fields), '<setup> : <sweep>_Incident'
    (direct Tx->Rx coupling only) and '<setup> : <sweep>_Scattered' (target only, the default,
    `scattered_solution`); target signatures need the scattered one.

    S is read as AEDT's real and imaginary parts and stored complex: that is AEDT's own form,
    exact, and free of the -inf dB of a zero channel and of phase wrapping (magnitude/phase data
    from elsewhere rebuild with `polarimetry_processor.complex_from_db_deg`).

    The .npz holds `s` (n_variations, n_expressions, n_frequencies) complex, `frequencies_hz`,
    `sweep_values`, `expressions` (AEDT names, e.g. 'S(HornCluster1_RxH,HornCluster1_TxV)', i.e.
    S(receive, transmit)), and `metadata` (JSON: project, design, setup, sweep variable and the
    nominal design variables; for an `sbr_simulator` scene also its boresight and the target's
    steering angles per variation, `scene_target`). Without a sweep variable only the nominal
    variation is written.
    When the linked antenna was built by `hfss_array_builder`, its port table (`ports_path`, or
    found next to the project) is added to the metadata as 'ports': per excitation the builder
    port, role, polarization, field sign and phase centre (mm, at build time), which
    `MimoPolarimetryProcessor` uses to correct the signs and as virtual-channel positions.
    The result data must exist: an SBR+ design whose results were invalidated (e.g. by an edit
    of the linked source) reports "No Data Available" until it is re-solved. Every variation must
    be solved on the same frequency grid; otherwise nothing is written
    (`polarimetry_processor.check_complete_sweep`). Variations with identical data are written
    with a warning (`polarimetry_processor.identical_variations`).
    """
    output_path = resolve_output_path(output_path)
    setup_sweep = setup_sweep or scattered_solution(app)
    print(f"Reading '{setup_sweep}' of '{app.design_name}'...")
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
    # AEDT returns all variations on the union of their frequency grids, with NaN where a
    # variation has no solution; such an export is refused rather than written
    check_complete_sweep(s, frequencies * scale, sweep_values, sweep_variable)
    if warning := identical_variations(s, sweep_values, sweep_variable):
        print(f"  Warning: {warning}")
    metadata = {
        "project": app.project_name,
        "design": app.design_name,
        "setup_sweep": setup_sweep,
        "sweep_variable": sweep_variable,
        "nominal_variables": app.available_variations.nominal_values,
        **scene_target(app, sweep_variable, sweep_values),
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
        default=None,
        help=(
            "Setup and sweep name (default: the design's '<setup> : <sweep>_Scattered', the "
            "target only)"
        ),
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

    # Fail on a bad output path before attaching to AEDT and reading the whole solution
    try:
        output_path = resolve_output_path(args.output)
    except (FileNotFoundError, PermissionError) as e:
        sys.exit(f"[ERROR] {e}")
    print(f"Exporting to {output_path}")
    try:
        desktop = connect_live_session()
    except (NativeAedtRefusedError, NoLiveSessionError) as e:
        sys.exit(f"[ERROR] {e}")
    try:
        require_designs(desktop, args.project, [args.design])
        app = Hfss(project=args.project, design=args.design, new_desktop=False, close_on_exit=False)
        export_sbr_parametric_sweep(
            app, output_path, args.setup, args.sweep_variable, ports_path=args.port_table
        )
    except IncompleteSweepError as e:
        sys.exit(f"[ERROR] Nothing written. {e}")
    except MissingDesignError as e:
        sys.exit(f"[ERROR] {e}")
    finally:
        desktop.release_desktop(close_projects=False, close_on_exit=False)


if __name__ == "__main__":
    main()
