#!/usr/bin/env python3
"""
Bistaticity study: polarimetric signatures of one target against a swept variable (e.g. the
Tx-Rx spacing), from an SBR+ sweep exported by `sbr_results` (scattered solution).

    bistaticity_analysis data/interim/bistatic_sweep/<sweep>.npz \
        [--out <folder>] [--azimuth DEG] [--elevation DEG] [--quad NAME]

All processing is `CanonicalTargetProcessor`'s: the target's scattering matrices beamformed
towards (azimuth, elevation) (boresight by default, where no geometry is needed), range profiles
and, for arrays with several virtual channels per polarization, azimuth-elevation maps. This
script only combines and plots them, and shows only target-agnostic quantities: no reference
target, no thresholds, no classification. The signatures are plotted twice: against the swept
variable (a line at the centre frequency and a min-max band over the sweep) and as frequency x
swept-variable maps, which show where in the band the spread lies. Figures and a CSV land in
<out>/<sweep>/, one folder per export (named after the .npz file, created if needed), with <out>
defaulting to the project's data/processed/bistaticity; the README ("Bistaticity study")
explains the signatures.
"""

import argparse
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.config import PROCESSED_DATA_DIR  # noqa: E402
from pomaar.simulator.array_templates import polarimetric_quads  # noqa: E402
from pomaar.simulator.polarimetry_processor import (  # noqa: E402
    POLARIMETRIC_CHANNELS,
    CanonicalTargetProcessor,
    IncompleteSweepError,
    bistatic_angle_deg,
    channel_ports,
    channel_powers,
    db10,
    identical_variations,
    pauli_powers,
    total_power,
)

DEFAULT_OUTPUT_ROOT = PROCESSED_DATA_DIR / "bistaticity"
SERIES_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e4e4e4"
PAULI_NAMES = ["k₁ ∝ HH+VV", "k₂ ∝ HH−VV", "k₃ ∝ HV+VH", "k₄ ∝ j(HV−VH)"]
# Phase differences (deg) are circular: their bands are taken around the centre-frequency value
PHASE_KEYS = {
    "copol_phase_difference_deg",
    "co_cross_phase_txh_deg",
    "co_cross_phase_txv_deg",
    "reciprocity_phase_deg",
}
# Bounded quantities are shown over their whole domain, so that noise is not magnified
DOMAINS = {"entropy": (0.0, 1.0), "anisotropy": (0.0, 1.0), "mean_alpha_deg": (0.0, 90.0)}


def style_axis(ax, title, ylabel):
    ax.set_title(title, loc="left", fontsize=9.5, color=INK)
    ax.set_ylabel(ylabel, fontsize=8.5, color=MUTED)
    ax.tick_params(labelsize=8, colors=MUTED)
    ax.grid(True, color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def wrap_deg(angle):
    """Angles (deg) wrapped to [-180, 180)."""
    return (np.asarray(angle) + 180.0) % 360.0 - 180.0


def band_limits(values, centre_index, circular=False):
    """
    Min-max over the band of `values` (n_var x n_freq). For a phase (deg, `circular`) it is the
    min-max of the deviation from the centre-frequency value, wrapped to +-180 deg and added to
    that value, so a phase scattered around +-180 deg does not span the whole circle.
    """
    if not circular:
        return values.min(axis=1), values.max(axis=1)
    line = values[:, centre_index]
    deviation = wrap_deg(values - line[:, None])
    return line + deviation.min(axis=1), line + deviation.max(axis=1)


def set_domain(ax, low, high, margin=0.03):
    pad = margin * (high - low)
    ax.set_ylim(low - pad, high + pad)


def band_plot(ax, x, values, centre_index, colour, label, circular=False):
    """
    Centre-frequency line plus min-max envelope over the band (values: n_var x n_freq). A
    circular band is drawn on a [-180, 180] deg axis and repeated 360 deg up and down, so the
    part that crosses +-180 deg reappears at the other edge.
    """
    low, high = band_limits(values, centre_index, circular)
    for shift in (-360.0, 0.0, 360.0) if circular else (0.0,):
        ax.fill_between(x, low + shift, high + shift, color=colour, alpha=0.18, lw=0)
    ax.plot(x, values[:, centre_index], color=colour, lw=2, marker="o", markersize=5, label=label)
    if circular:
        ax.set_ylim(-180, 180)
        ax.set_yticks(np.arange(-180, 181, 90))


def plot_range_profiles(out_dir, maps, x, x_label, picks):
    """Beamformed range profiles per polarimetric channel, relative to the total-power peak."""
    powers = channel_powers(maps["channels"])
    total = total_power(maps["channels"])
    reference = total.reshape(total.shape[0], -1).max(axis=-1)[:, None, None]
    fig, axes = plt.subplots(1, 5, figsize=(22, 4.0), constrained_layout=True, sharey=True)
    for ax, name in zip(axes, ["total"] + list(POLARIMETRIC_CHANNELS)):
        values = total / reference if name == "total" else powers[name]
        for colour, index in zip(SERIES_COLOURS, picks):
            ax.plot(
                maps["range_m"],
                db10(values[index, :, 0]),
                color=colour,
                lw=1.4,
                label=f"{x_label} = {x[index]:g}",
            )
        style_axis(ax, "Total power" if name == "total" else f"|b_{name}|²", "dB re total peak")
        ax.set_xlabel("range (m)", fontsize=8.5, color=MUTED)
        ax.set_ylim(-100, 5)
    axes[0].legend(fontsize=7.5, frameon=False)
    fig.suptitle(
        f"Range profiles towards the target ({maps['azimuth_deg'][0]:g}°, "
        f"{maps['elevation_deg']:g}°); unambiguous range {maps['range_m'][-1]:.2f} m",
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(os.path.join(out_dir, "range_profiles.png"), dpi=150)
    plt.close(fig)


def plot_maps(out_dir, maps, index, x_label, x_value):
    """
    Azimuth-elevation maps of the total power and of every Pauli power |k_i|^2 at the target
    range, all relative to the peak span (so that each map shows where its power is).
    """
    total = total_power(maps["channels"])[index]
    pauli = pauli_powers(maps["channels"])[index]
    peak_span = pauli.sum(axis=-1).max()
    extent = [
        maps["azimuth_deg"][0],
        maps["azimuth_deg"][-1],
        maps["elevation_deg"][0],
        maps["elevation_deg"][-1],
    ]
    fig, axes = plt.subplots(1, 5, figsize=(22, 4.0), constrained_layout=True)
    images = [(db10(total / total.max()), "Total power (dB re peak)", (-40, 0))] + [
        (db10(pauli[..., i] / peak_span), f"|{name.split()[0]}|² (dB re peak span)", (-40, 0))
        for i, name in enumerate(PAULI_NAMES)
    ]
    for ax, (image, title, limits) in zip(axes, images):
        shown = ax.imshow(image.T, origin="lower", extent=extent, aspect="auto", cmap="viridis")
        shown.set_clim(*limits)
        fig.colorbar(shown, ax=ax, shrink=0.85)
        style_axis(ax, title, "elevation (deg)")
        ax.set_xlabel("azimuth (deg)", fontsize=8.5, color=MUTED)
    fig.suptitle(
        f"Azimuth-elevation maps at {maps['range_m'][index]:.3f} m, {x_label} = {x_value:g}",
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(os.path.join(out_dir, "azimuth_elevation_maps.png"), dpi=150)
    plt.close(fig)


def plot_signatures(
    out_dir, features, x, x_label, centre_index, frequencies, title, name="signatures_vs_sweep"
):
    """
    Target-agnostic signatures against the swept variable: a line at the centre frequency and,
    for per-frequency quantities, the min-max band over the sweep (around the line for phases).
    Phases span [-180, 180] deg, Pauli fractions end at 0 dB and H, A and mean alpha span their
    whole domain; the other axes fit the data.
    """
    fig, axes = plt.subplots(3, 4, figsize=(20, 11.5), constrained_layout=True)
    axes = axes.ravel()
    panels = [
        (
            "(1) Channel returns |S_pq|²",
            "dB",
            [(f"power_{q.lower()}_db", f"|S_{q}|²") for q in POLARIMETRIC_CHANNELS],
        ),
        ("(2) Span ‖S‖²", "dB", [("span_db", None)]),
        ("(3) Co-pol ratio |S_VV / S_HH|", "dB", [("copol_ratio_db", None)]),
        (
            "(4) Co-pol phase difference ∠(S_VV S_HH*)",
            "deg",
            [("copol_phase_difference_deg", None)],
        ),
        (
            "(5) Co- to cross-pol ratio per Tx polarization",
            "dB",
            [
                ("co_minus_cross_txh_db", "Tx H: |S_HH / S_VH|"),
                ("co_minus_cross_txv_db", "Tx V: |S_VV / S_HV|"),
            ],
        ),
        (
            "(6) Co- to cross-pol phase difference per Tx polarization",
            "deg",
            [
                ("co_cross_phase_txh_deg", "Tx H: ∠(S_HH S_VH*)"),
                ("co_cross_phase_txv_deg", "Tx V: ∠(S_VV S_HV*)"),
            ],
        ),
        ("(7) Cross-pol ratio |S_HV / S_VH|", "dB", [("reciprocity_db", None)]),
        ("(8) Cross-pol phase difference ∠(S_HV S_VH*)", "deg", [("reciprocity_phase_deg", None)]),
        ("(9) Normalized asymmetry |S_HV − S_VH| / ‖S‖", "dB", [("asymmetry_db", None)]),
        (
            "(10) Bistatic Pauli power fractions |k_i|² / ‖k‖²",
            "dB",
            [(("pauli_fractions_db", i), name) for i, name in enumerate(PAULI_NAMES)],
        ),
    ]
    for ax, (panel_title, unit, series) in zip(axes, panels):
        for colour, (key, label) in zip(SERIES_COLOURS, series):
            values = features[key[0]][..., key[1]] if isinstance(key, tuple) else features[key]
            band_plot(ax, x, values, centre_index, colour, label, circular=key in PHASE_KEYS)
        style_axis(ax, panel_title, unit)
        if series[0][0][0] == "pauli_fractions_db":
            low = ax.get_ylim()[0]
            set_domain(ax, low, 0.0)
        if len(series) > 1:
            ax.legend(fontsize=7, frameon=False)

    # H and A share the domain [0, 1], so they share a panel
    cloude = [
        (
            "(11) Cloude–Pottier entropy H and anisotropy A (4×4, band-averaged)",
            "",
            [("entropy", "H"), ("anisotropy", "A")],
        ),
        ("(12) Cloude–Pottier mean α (4×4, band-averaged)", "deg", [("mean_alpha_deg", None)]),
    ]
    for ax, (panel_title, unit, series) in zip(axes[10:], cloude):
        for colour, (key, label) in zip(SERIES_COLOURS, series):
            ax.plot(x, features[key], color=colour, lw=2, marker="o", markersize=5, label=label)
        style_axis(ax, panel_title, unit)
        set_domain(ax, *DOMAINS[series[0][0]])
        if len(series) > 1:
            ax.legend(fontsize=7, frameon=False)

    for ax in axes:
        ax.set_xlabel(x_label, fontsize=8.5, color=MUTED)
    fig.suptitle(
        f"{title}. Line: {frequencies[centre_index] / 1e9:.2f} GHz; band: min–max over "
        f"{frequencies[0] / 1e9:.2f}–{frequencies[-1] / 1e9:.2f} GHz (phases: of the wrapped "
        "deviation from the line).",
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(os.path.join(out_dir, f"{name}.png"), dpi=150)
    plt.close(fig)


def cell_edges(centres):
    """Cell edges around (sorted) cell centres, for pcolormesh."""
    centres = np.asarray(centres, dtype=float)
    if centres.size == 1:
        return np.array([centres[0] - 0.5, centres[0] + 0.5])
    middle = (centres[1:] + centres[:-1]) / 2
    return np.concatenate([[2 * centres[0] - middle[0]], middle, [2 * centres[-1] - middle[-1]]])


def finite_limits(*arrays):
    finite = np.concatenate([np.asarray(a)[np.isfinite(a)] for a in arrays])
    return (finite.min(), finite.max()) if finite.size else (-1.0, 1.0)


def plot_signature_maps(
    out_dir, features, x, x_label, frequencies, title, name="signatures_vs_frequency"
):
    """
    Every per-frequency signature as a frequency x swept-variable map: the unsummarized data
    behind the bands of `plot_signatures`. Every panel has its own colour scale over its data's
    min and max, so the variation within each quantity stays visible; ratios use a scale
    symmetric about 0 dB. Phases are shown as their circular mean plus the wrapped deviation
    from it, on a cyclic scale when the deviations reach around the circle.
    """
    order = np.argsort(x)
    x_edges = cell_edges(np.asarray(x)[order])
    f_edges = cell_edges(np.asarray(frequencies) / 1e9)
    channels = [(f"power_{q.lower()}_db", f"|S_{q}|² (dB)", None) for q in POLARIMETRIC_CHANNELS]
    panels = (
        channels
        + [
            ("span_db", "Span ‖S‖² (dB)", None),
            ("copol_ratio_db", "Co-pol ratio |S_VV / S_HH| (dB)", "ratio"),
            ("copol_phase_difference_deg", "Co-pol phase difference ∠(S_VV S_HH*) (deg)", "phase"),
            ("co_minus_cross_txh_db", "Tx H: |S_HH / S_VH| (dB)", None),
            ("co_cross_phase_txh_deg", "Tx H: ∠(S_HH S_VH*) (deg)", "phase"),
            ("co_minus_cross_txv_db", "Tx V: |S_VV / S_HV| (dB)", None),
            ("co_cross_phase_txv_deg", "Tx V: ∠(S_VV S_HV*) (deg)", "phase"),
            ("reciprocity_db", "Cross-pol ratio |S_HV / S_VH| (dB)", "ratio"),
            ("reciprocity_phase_deg", "Cross-pol phase difference ∠(S_HV S_VH*) (deg)", "phase"),
            ("asymmetry_db", "Normalized asymmetry |S_HV − S_VH| / ‖S‖ (dB)", None),
        ]
        + [
            (("pauli_fractions_db", i), f"Pauli fraction |{name.split()[0]}|² / ‖k‖² (dB)", None)
            for i, name in enumerate(PAULI_NAMES)
        ]
    )

    def values_of(key):
        values = features[key[0]][..., key[1]] if isinstance(key, tuple) else features[key]
        return values[order]

    n_columns = 4
    n_rows = -(-len(panels) // n_columns)
    fig, axes = plt.subplots(
        n_rows, n_columns, figsize=(20, 3.6 * n_rows), constrained_layout=True, sharey=True
    )
    for ax, (key, panel_title, group) in zip(axes.ravel(), panels):
        values = values_of(key)
        if group == "phase":
            mean = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(values)))))
            values = mean + wrap_deg(values - mean)
            extent = np.max(np.abs(values - mean))
            if extent > 170.0:
                cmap, limits = "twilight", (mean - 180.0, mean + 180.0)
            else:
                cmap, limits = "viridis", (mean - extent, mean + extent)
        elif group == "ratio":
            extent = np.max(np.abs(finite_limits(values)))
            cmap, limits = "RdBu_r", (-extent, extent)
        else:
            cmap, limits = "viridis", finite_limits(values)
        shown = ax.pcolormesh(
            x_edges, f_edges, values.T, cmap=cmap, vmin=limits[0], vmax=limits[1], rasterized=True
        )
        fig.colorbar(shown, ax=ax, shrink=0.9)
        style_axis(ax, panel_title, "frequency (GHz)")
        ax.grid(False)
        ax.set_xlabel(x_label, fontsize=8.5, color=MUTED)
    for ax in axes.ravel()[len(panels) :]:
        ax.set_visible(False)
    fig.suptitle(
        f"{title}. Signatures over frequency and {x_label}.",
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(os.path.join(out_dir, f"{name}.png"), dpi=150)
    plt.close(fig)


def write_csv(out_dir, features, x, x_label, centre_index, name="signatures_vs_sweep", extra=None):
    """
    Every signature: per-frequency ones at the centre frequency and their band min/max (for
    phases, the line plus the min/max of the wrapped deviation from it, as plotted). `extra`
    ({column: values}) adds columns after the first, e.g. where each row comes from.
    """
    extra = extra or {}
    per_frequency = [k for k, v in features.items() if np.ndim(v) == 2]
    per_variation = [k for k, v in features.items() if np.ndim(v) == 1]
    bands = {k: band_limits(features[k], centre_index, k in PHASE_KEYS) for k in per_frequency}
    with open(os.path.join(out_dir, f"{name}.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [x_label]
            + list(extra)
            + [f"{k}@f0" for k in per_frequency]
            + [f"{k}_band_{e}" for k in per_frequency for e in ("min", "max")]
            + [f"pauli_k{i}_fraction_db@f0" for i in range(1, 5)]
            + per_variation
        )
        for index, value in enumerate(x):
            writer.writerow(
                [value]
                + [column[index] for column in extra.values()]
                + [features[k][index, centre_index] for k in per_frequency]
                + [bands[k][e][index] for k in per_frequency for e in (0, 1)]
                + list(features["pauli_fractions_db"][index, centre_index])
                + [features[k][index] for k in per_variation]
            )


def receiver_bistatic_angles(metadata, quad_ports):
    """
    Bistatic angle (deg) of one polarimetric quad per variation, from the scene geometry in the
    export's metadata (`positions_m` per port and `target_m`, per variation, in metres): the
    angle at the target between the mean Tx and the mean Rx position of the quad. None when the
    export carries no geometry.
    """
    positions, target = metadata.get("positions_m"), metadata.get("target_m")
    if not positions or not target:
        return None
    tx = np.mean([positions[quad_ports["TxH"]], positions[quad_ports["TxV"]]], axis=0)
    rx = np.mean([positions[quad_ports["RxH"]], positions[quad_ports["RxV"]]], axis=0)
    return bistatic_angle_deg(tx, rx, target)


def pool_receivers(path, azimuth=0.0, elevation=0.0, tolerance_deg=1e-6):
    """
    The signatures of every receiver (polarimetric quad) and variation of an export, pooled and
    sorted by bistatic angle: (features, angles, receivers, sweep values). Points of equal angle
    (within `tolerance_deg`, e.g. the monostatic receiver at every variation) are kept once; the
    largest difference of their per-frequency signatures is reported as a consistency check.
    """
    expressions = np.load(path)["expressions"]
    quads = polarimetric_quads([port for pair in channel_ports(expressions) for port in pair])
    rows, angles, receivers, sweep_values = [], [], [], []
    for quad, quad_ports in quads.items():
        processor = CanonicalTargetProcessor.from_npz(path, quad=quad)
        beta = receiver_bistatic_angles(processor.metadata, quad_ports)
        if beta is None:
            sys.exit(
                f"[ERROR] {os.path.basename(path)} has no scene geometry (metadata positions_m "
                "and target_m), so receivers cannot be placed by bistatic angle."
            )
        rows.append(processor.signatures(matrices=processor.target_matrices(azimuth, elevation)))
        angles.append(beta)
        receivers += [quad.split("/")[0]] * len(beta)
        sweep_values.append(processor.sweep_values)
    features = {k: np.concatenate([r[k] for r in rows]) for k in rows[0]}
    angles, sweep_values = np.concatenate(angles), np.concatenate(sweep_values)
    order = np.argsort(angles, kind="stable")
    keep, spread = [], 0.0
    per_frequency = [k for k, v in features.items() if np.ndim(v) == 2]
    for index in order:
        if keep and abs(angles[index] - angles[keep[-1]]) < tolerance_deg:
            for k in per_frequency:
                difference = np.abs(features[k][index] - features[k][keep[-1]])
                if "phase" in k:
                    difference = np.abs(wrap_deg(difference))
                spread = max(spread, float(np.nanmax(difference)))
            continue
        keep.append(index)
    print(
        f"{len(quads)} receivers x {len(order) // len(quads)} variations -> {len(keep)} distinct "
        f"bistatic angles ({angles[keep].min():g} ... {angles[keep].max():g} deg); largest "
        f"difference between points of equal angle: {spread:.3g} (dB or deg)"
    )
    keep = np.asarray(keep)
    return (
        {k: v[keep] for k, v in features.items()},
        angles[keep],
        [receivers[i] for i in keep],
        sweep_values[keep],
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("sweep", help=".npz written by pomaar.simulator.sbr_results")
    parser.add_argument(
        "--out",
        default=str(DEFAULT_OUTPUT_ROOT),
        help=(
            "Parent folder; the results go to <out>/<export file name>/ "
            f"(default: {DEFAULT_OUTPUT_ROOT})"
        ),
    )
    parser.add_argument(
        "--azimuth", type=float, default=0.0, help="Steering azimuth (deg, towards +x)"
    )
    parser.add_argument(
        "--elevation", type=float, default=0.0, help="Steering elevation (deg, towards +y)"
    )
    parser.add_argument(
        "--quad", default=None, help="Keep one polarimetric quad (e.g. 'Rx_1/Tx'); default all"
    )
    parser.add_argument(
        "--window", default="hann", help="Frequency window of the range profiles (default: hann)"
    )
    parser.add_argument(
        "--per-receiver",
        action="store_true",
        help=(
            "Treat every receiver (polarimetric quad) separately and plot all of them, pooled, "
            "against their bistatic angle (needs the scene geometry in the export's metadata)"
        ),
    )
    args = parser.parse_args(argv)

    # Check the input and the output folder before any processing
    sweep = os.path.abspath(os.path.expanduser(args.sweep))
    # One folder per export, so that analysing another data set never overwrites results
    stem = os.path.splitext(os.path.basename(sweep))[0]
    out = os.path.join(os.path.abspath(os.path.expanduser(args.out)), stem)
    if not os.path.isfile(sweep):
        sys.exit(f"[ERROR] No export at {sweep} (relative paths resolve against {os.getcwd()}).")
    try:
        os.makedirs(out, exist_ok=True)
    except OSError as e:
        sys.exit(f"[ERROR] Cannot create the output folder {out} ({e}).")
    if not os.access(out, os.W_OK):
        sys.exit(f"[ERROR] The output folder {out} is not writable.")
    args.sweep, args.out = sweep, out

    try:
        processor = CanonicalTargetProcessor.from_npz(args.sweep, quad=args.quad)
    except IncompleteSweepError as e:
        sys.exit(f"[ERROR] {os.path.basename(sweep)}: {e}")
    frequencies = processor.frequencies_hz
    x_label = processor.metadata.get("sweep_variable") or "variation"
    x = processor.sweep_values
    if warning := identical_variations(processor.raw()["s"], x, x_label):
        print(f"Warning: {warning}")
    if np.all(np.isnan(x)):
        x = np.arange(len(x), dtype=float)
    centre_index = len(frequencies) // 2
    design = processor.metadata.get("design") or os.path.basename(args.sweep)
    print(
        f"{design} ({processor.metadata.get('setup_sweep')}): {len(x)} values of {x_label}, "
        f"{len(frequencies)} points, {frequencies[0] / 1e9:.2f}-{frequencies[-1] / 1e9:.2f} GHz; "
        "virtual channels "
        + ", ".join(f"{q} {len(i)}" for q, i in processor.polarimetric_channels.items())
    )
    off_boresight = args.azimuth != 0.0 or args.elevation != 0.0
    if off_boresight and not processor.positions_known:
        print(
            "Warning: the export has no port positions (builder port table), so the channels "
            "are not steered off boresight."
        )

    if args.per_receiver:
        features, angles, receivers, values = pool_receivers(
            args.sweep, args.azimuth, args.elevation
        )
        title = (
            f"{design}, every receiver at its bistatic angle, steered to ({args.azimuth:g}°, "
            f"{args.elevation:g}°)"
        )
        label = "bistatic angle (deg)"
        name = "signatures_vs_bistatic_angle"
        plot_signatures(args.out, features, angles, label, centre_index, frequencies, title, name)
        plot_signature_maps(
            args.out, features, angles, label, frequencies, title, f"{name}_and_frequency"
        )
        write_csv(
            args.out,
            features,
            angles,
            label,
            centre_index,
            name,
            extra={"receiver": receivers, x_label: values},
        )
        print(f"Figures and CSV written to {args.out}")
        return

    matrices = processor.target_matrices(args.azimuth, args.elevation)
    features = processor.signatures(matrices=matrices)

    picks = sorted(set(np.linspace(0, len(x) - 1, min(4, len(x))).round().astype(int)))
    profiles = processor.range_azimuth_maps(
        [args.azimuth], elevation_deg=args.elevation, window=args.window
    )
    plot_range_profiles(args.out, profiles, x, x_label, picks)
    if max(len(i) for i in processor.polarimetric_channels.values()) > 1:
        maps = processor.azimuth_elevation_maps(
            np.arange(-60, 60.5, 1.0), np.arange(-30, 30.5, 1.0), window=args.window
        )
        plot_maps(args.out, maps, 0, x_label, x[0])
    else:
        print("One virtual channel per polarization: no angular resolution, so no maps.")
    title = f"{design}, steered to ({args.azimuth:g}°, {args.elevation:g}°)"
    plot_signatures(args.out, features, x, x_label, centre_index, frequencies, title)
    plot_signature_maps(args.out, features, x, x_label, frequencies, title)
    write_csv(args.out, features, x, x_label, centre_index)
    print(f"Figures and CSV written to {args.out}")


if __name__ == "__main__":
    main()
