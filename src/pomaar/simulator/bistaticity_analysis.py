#!/usr/bin/env python3
"""
Bistaticity study: polarimetric signatures of one target against a swept variable (e.g. the
Tx-Rx spacing), from an SBR+ sweep exported by `sbr_results` (scattered solution).

    bistaticity_analysis data/interim/bistatic_sweep/<sweep>.npz \
        [--out data/processed/bistaticity] [--azimuth DEG] [--elevation DEG] [--quad NAME]

All processing is `CanonicalTargetProcessor`'s: the target's scattering matrices beamformed
towards (azimuth, elevation) (boresight by default, where no geometry is needed), range profiles
and, for arrays with several virtual channels per polarization, azimuth-elevation maps. This
script only combines and plots them, and shows only target-agnostic quantities: no reference
target, no thresholds, no classification. Figures and a CSV land in --out; the README
("Bistaticity study") explains the signatures.
"""

import argparse
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.simulator.polarimetry_processor import (  # noqa: E402
    POLARIMETRIC_CHANNELS,
    CanonicalTargetProcessor,
    channel_powers,
    db10,
    pauli_powers,
    total_power,
)

SERIES_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e4e4e4"
PAULI_NAMES = ["k₁ ∝ HH+VV", "k₂ ∝ HH−VV", "k₃ ∝ HV+VH", "k₄ ∝ j(HV−VH)"]


def style_axis(ax, title, ylabel):
    ax.set_title(title, loc="left", fontsize=9.5, color=INK)
    ax.set_ylabel(ylabel, fontsize=8.5, color=MUTED)
    ax.tick_params(labelsize=8, colors=MUTED)
    ax.grid(True, color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def band_plot(ax, x, values, centre_index, colour, label):
    """Centre-frequency line plus min-max envelope over the band (values: n_var x n_freq)."""
    ax.fill_between(x, values.min(axis=1), values.max(axis=1), color=colour, alpha=0.18, lw=0)
    ax.plot(x, values[:, centre_index], color=colour, lw=2, marker="o", markersize=5, label=label)


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


def plot_signatures(out_dir, features, x, x_label, centre_index, frequencies, title):
    """
    Target-agnostic signatures against the swept variable: a line at the centre frequency and,
    for per-frequency quantities, the min-max band over the sweep.
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
            "(5) Co-pol minus cross-pol return per Tx polarization",
            "dB",
            [
                ("co_minus_cross_txh_db", "Tx H: |S_HH| − |S_VH|"),
                ("co_minus_cross_txv_db", "Tx V: |S_VV| − |S_HV|"),
            ],
        ),
        ("(6) Cross-pol ratio |S_HV / S_VH|", "dB", [("reciprocity_db", None)]),
        ("(7) Cross-pol phase difference ∠(S_HV S_VH*)", "deg", [("reciprocity_phase_deg", None)]),
        ("(8) Normalized asymmetry |S_HV − S_VH| / ‖S‖", "dB", [("asymmetry_db", None)]),
        (
            "(9) Bistatic Pauli power fractions |k_i|² / ‖k‖²",
            "dB",
            [(("pauli_fractions_db", i), name) for i, name in enumerate(PAULI_NAMES)],
        ),
    ]
    for ax, (panel_title, unit, series) in zip(axes, panels):
        for colour, (key, label) in zip(SERIES_COLOURS, series):
            values = features[key[0]][..., key[1]] if isinstance(key, tuple) else features[key]
            band_plot(ax, x, values, centre_index, colour, label)
        style_axis(ax, panel_title, unit)
        if len(series) > 1:
            ax.legend(fontsize=7, frameon=False)

    cloude = [
        ("(10) Cloude–Pottier entropy H (4×4, band-averaged)", "H", "entropy"),
        ("(11) Cloude–Pottier anisotropy A (4×4, band-averaged)", "A", "anisotropy"),
        ("(12) Cloude–Pottier mean α (4×4, band-averaged)", "deg", "mean_alpha_deg"),
    ]
    for ax, (panel_title, unit, key) in zip(axes[9:], cloude):
        ax.plot(x, features[key], color=SERIES_COLOURS[0], lw=2, marker="o", markersize=5)
        style_axis(ax, panel_title, unit)

    for ax in axes:
        ax.set_xlabel(x_label, fontsize=8.5, color=MUTED)
    fig.suptitle(
        f"{title}. Line: {frequencies[centre_index] / 1e9:.2f} GHz; band: min–max over "
        f"{frequencies[0] / 1e9:.2f}–{frequencies[-1] / 1e9:.2f} GHz.",
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(os.path.join(out_dir, "signatures_vs_sweep.png"), dpi=150)
    plt.close(fig)


def write_csv(out_dir, features, x, x_label, centre_index):
    """Every signature: per-frequency ones at the centre frequency and their band min/max."""
    per_frequency = [k for k, v in features.items() if np.ndim(v) == 2]
    per_variation = [k for k, v in features.items() if np.ndim(v) == 1]
    with open(os.path.join(out_dir, "signatures_vs_sweep.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [x_label]
            + [f"{k}@f0" for k in per_frequency]
            + [f"{k}_band_{e}" for k in per_frequency for e in ("min", "max")]
            + [f"pauli_k{i}_fraction_db@f0" for i in range(1, 5)]
            + per_variation
        )
        for index, value in enumerate(x):
            writer.writerow(
                [value]
                + [features[k][index, centre_index] for k in per_frequency]
                + [f(features[k][index]) for k in per_frequency for f in (np.min, np.max)]
                + list(features["pauli_fractions_db"][index, centre_index])
                + [features[k][index] for k in per_variation]
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("sweep", help=".npz written by pomaar.simulator.sbr_results")
    parser.add_argument("--out", default="data/processed/bistaticity")
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
    args = parser.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)

    processor = CanonicalTargetProcessor.from_npz(args.sweep, quad=args.quad)
    frequencies = processor.frequencies_hz
    x_label = processor.metadata.get("sweep_variable") or "variation"
    x = processor.sweep_values
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
    write_csv(args.out, features, x, x_label, centre_index)
    print(f"Figures and CSV written to {args.out}")


if __name__ == "__main__":
    main()
