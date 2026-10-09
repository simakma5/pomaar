#!/usr/bin/env python3
"""
Bistaticity study of the dual-polarised horn cluster: polarimetric signatures versus the
co-polar Tx-Rx spacing, from an SBR+ sweep exported by `sbr_results`.

    bistaticity_analysis data/interim/bistatic_sweep/<sweep>.npz \
        [--out data/processed/bistaticity] [--gate auto|on|off] [--background <no_target>.npz]

Figures and a metrics CSV land in --out. The README ("Bistaticity study") explains every step.
"""

import argparse
import copy
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.simulator.polarimetry_processor import (  # noqa: E402
    REFERENCE_SPHERE,
    SPEED_OF_LIGHT,
    MimoPolarimetryProcessor,
    bistatic_path,
    delay_profile,
    delay_window_peak,
    pauli_power_fractions,
    plane_wave_phase,
    spherical_wave_phase,
    wrap_path,
)

GATE_WINDOW = ("kaiser", 14)
# Guard band beside the gate (offsets from the target path, m), on the side away from the direct
# coupling, and the minimum target-to-guard ratio for a co-polar point to count as clean
GUARD_OFFSETS = (-0.40, -0.18)
MIN_GATE_SNR_DB = 20.0
# A target peak displaced from the sweep's median target path by more than this is contaminated
MAX_PATH_DEVIATION_M = 0.005
# Cross-pol leakage counts as resolved only this far above its measurement floor
MIN_CROSSPOL_MARGIN_DB = 10.0
SERIES_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e4e4e4"


def cluster_positions(copolar_spacing_m):
    """
    Phase-centre positions (m) of the DualPolHornCluster elements, mirroring its CS origins.

    The four elements sit on a square of diagonal `copolar_spacing_m` around the origin: co-polar
    Tx and Rx on opposite corners, so TxH-RxH and TxV-RxV are one spacing apart along the two
    diagonals and every cross-polar pair is spacing / sqrt(2) apart along x.
    """
    offset = copolar_spacing_m / (2 * np.sqrt(2))
    return {
        "TxV": np.array([offset, offset, 0.0]),
        "RxV": np.array([-offset, -offset, 0.0]),
        "TxH": np.array([offset, -offset, 0.0]),
        "RxH": np.array([-offset, offset, 0.0]),
    }


def pair_positions(positions, rx, tx):
    """Returns (tx_position, rx_position) of the channel S[rx, tx]."""
    names = "HV"
    return positions[f"Tx{names[tx]}"], positions[f"Rx{names[rx]}"]


def locate_target_path(matrices, frequencies, expected_path, search_half_width):
    """
    Finds the two-way path of the target peak, common to all four channels.

    `matrices` is (n_freq, 2, 2) for one variation. Searches the zero-padded delay profile of
    |S_HH|^2 + |S_VV|^2 within +-search_half_width of the (wrapped) expected path, so the
    co-polar returns set one gate for every channel.
    """
    copolar = np.moveaxis(matrices[..., [0, 1], [0, 1]], 0, -1)  # (2, n_freq)
    path, profile = delay_profile(
        copolar, frequencies, window=GATE_WINDOW, n_fft=16 * len(frequencies)
    )
    power = np.sum(np.abs(profile) ** 2, axis=0)
    unambiguous = SPEED_OF_LIGHT / (frequencies[1] - frequencies[0])
    expected_wrapped = wrap_path(expected_path, frequencies)
    distance = np.abs(wrap_path(path - expected_wrapped, frequencies))
    candidates = np.flatnonzero(distance <= search_half_width)
    peak = candidates[np.argmax(power[candidates])]
    # Unwrap the found peak next to the expected absolute path
    return expected_path + wrap_path(path[peak] - expected_wrapped, frequencies), unambiguous


def db20(x):
    return 20 * np.log10(np.abs(x) + 1e-300)


def db10(x):
    return 10 * np.log10(np.abs(x) + 1e-300)


def style_axis(ax, title, ylabel):
    ax.set_title(title, loc="left", fontsize=9.5, color=INK)
    ax.set_ylabel(ylabel, fontsize=8.5, color=MUTED)
    ax.tick_params(labelsize=8, colors=MUTED)
    ax.grid(True, color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def band_plot(ax, x, values, valid, centre_index, colour, label, flagged=None):
    """
    Centre-frequency line plus min-max envelope over the valid band. Spacings in `flagged`
    (contaminated gate) get hollow markers.
    """
    band = values[:, valid]
    ax.fill_between(x, band.min(axis=1), band.max(axis=1), color=colour, alpha=0.18, linewidth=0)
    ax.plot(
        x, values[:, centre_index], color=colour, linewidth=2, marker="o", markersize=5, label=label
    )
    if flagged is not None and np.any(flagged):
        ax.plot(
            x[flagged],
            values[flagged, centre_index],
            linestyle="none",
            marker="o",
            markersize=5,
            markerfacecolor="white",
            markeredgecolor=colour,
            markeredgewidth=1.5,
        )


def add_bistatic_axis(ax, spacing_lambda, wavelength, target_distance):
    """Secondary x-axis with the co-polar bistatic angle beta = 2 atan(d / 2R)."""

    def to_beta(x):
        return np.degrees(2 * np.arctan(np.asarray(x) * wavelength / (2 * target_distance)))

    def to_spacing(beta):
        return 2 * target_distance * np.tan(np.radians(np.asarray(beta)) / 2) / wavelength

    top = ax.secondary_xaxis("top", functions=(to_beta, to_spacing))
    top.tick_params(labelsize=7, colors=MUTED)
    top.set_xlabel("bistatic angle β (deg)", fontsize=7.5, color=MUTED)


def plot_metrics(signatures):
    """Converts `MimoPolarimetryProcessor.signatures` to the dB / degree values plotted here."""
    return {
        "copol_ratio_db": db20(signatures["copol_ratio"]),
        "copol_phase_deg": np.degrees(np.angle(signatures["copol_ratio"])),
        "xpol_txh_db": db20(signatures["xpol_txh"]),
        "xpol_txv_db": db20(signatures["xpol_txv"]),
        "reciprocity_db": db20(signatures["reciprocity"]),
        "reciprocity_phase_deg": np.degrees(np.angle(signatures["reciprocity"])),
        "asymmetry_db": db20(signatures["asymmetry"]),
        "pauli_db": db10(signatures["pauli_fractions"]),
        "infidelity_sphere_db": db10(signatures["infidelity"]),
        "alpha_deg": signatures["alpha_deg"],
        "entropy_bistatic": signatures["entropy_bistatic"],
        "entropy_monostatic": signatures["entropy_monostatic"],
        "mean_alpha_bistatic_deg": signatures["mean_alpha_bistatic_deg"],
    }


def plot_validation(
    out_dir,
    s_matrices,
    gated,
    frequencies,
    spacing_lambda,
    target_paths,
    gate_half_width,
    valid,
    centre_index,
    gate_snr_db,
    flagged,
    gating,
):
    """Raw delay profiles with the gate, processed channel levels, and gate SNR versus spacing."""
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.2), constrained_layout=True)
    picks = [0, 2, len(spacing_lambda) // 2, len(spacing_lambda) - 1]
    for ax, (rx, tx, name) in zip(axes[:2], [(0, 0, "HH"), (0, 1, "HV")]):
        for colour, index in zip(SERIES_COLOURS, picks):
            path, profile = delay_profile(
                s_matrices[index, :, rx, tx],
                frequencies,
                window=GATE_WINDOW,
                n_fft=16 * len(frequencies),
            )
            ax.plot(
                path,
                db20(profile),
                color=colour,
                linewidth=1.4,
                label=f"d = {spacing_lambda[index]:.2f}λ",
            )
            centre = wrap_path(target_paths[index], frequencies)
            ax.axvspan(
                centre - gate_half_width,
                centre + gate_half_width,
                color=colour,
                alpha=0.08,
                linewidth=0,
            )
        ax.axvspan(
            centre + GUARD_OFFSETS[0],
            centre + GUARD_OFFSETS[1],
            color=MUTED,
            alpha=0.06,
            linewidth=0,
            hatch="//",
        )
        ax.set_xlim(-0.8, 0.6)
        style_axis(ax, f"Raw delay profile, S_{name} (Kaiser β=14)", "|profile| (dB)")
        ax.set_xlabel("two-way path, wrapped (m)", fontsize=8.5, color=MUTED)
        ax.legend(fontsize=7.5, frameon=False)
    ax = axes[2]
    for colour, (rx, tx, name) in zip(
        SERIES_COLOURS, [(0, 0, "HH"), (1, 1, "VV"), (0, 1, "HV"), (1, 0, "VH")]
    ):
        band_plot(
            ax,
            spacing_lambda,
            db20(gated[:, :, rx, tx]),
            valid,
            centre_index,
            colour,
            f"|S_{name}|",
            flagged,
        )
    style_axis(
        ax,
        "Gated target return (hollow: contaminated gate)" if gating else "Target return (ungated)",
        "|S| (dB)",
    )
    ax.set_xlabel("co-polar spacing d (λ₀)", fontsize=8.5, color=MUTED)
    ax.legend(fontsize=7.5, frameon=False)
    ax = axes[3]
    for colour, (rx, tx, name) in zip(
        SERIES_COLOURS, [(0, 0, "HH"), (1, 1, "VV"), (0, 1, "HV"), (1, 0, "VH")]
    ):
        ax.plot(
            spacing_lambda,
            gate_snr_db[:, rx, tx],
            color=colour,
            linewidth=2,
            marker="o",
            markersize=5,
            label=name,
        )
    ax.axhline(MIN_GATE_SNR_DB, color=MUTED, linewidth=0.8, linestyle="--")
    style_axis(ax, "Gate peak / guard-band peak" + ("" if gating else " (info only)"), "dB")
    ax.set_xlabel("co-polar spacing d (λ₀)", fontsize=8.5, color=MUTED)
    ax.legend(fontsize=7.5, frameon=False)
    fig.savefig(os.path.join(out_dir, "validation_gating.png"), dpi=150)
    plt.close(fig)


def plot_signatures(
    out_dir,
    metrics,
    spacing_lambda,
    valid,
    centre_index,
    wavelength,
    target_distance,
    frequencies,
    flagged,
    flagged_crosspol,
    crosspol_floor_db,
    gating,
):
    fig, axes = plt.subplots(3, 3, figsize=(15, 11.5), constrained_layout=True)
    axes = axes.ravel()
    x = spacing_lambda
    panels = [
        ("(1) Co-pol imbalance |S_VV / S_HH|", "dB", [("copol_ratio_db", None)]),
        ("(2) Co-pol phase difference ∠(S_VV S_HH*)", "deg", [("copol_phase_deg", None)]),
        (
            "(3) Cross-pol leakage per Tx polarisation",
            "dB",
            [("xpol_txh_db", "Tx H: |S_VH / S_HH|"), ("xpol_txv_db", "Tx V: |S_HV / S_VV|")],
        ),
        ("(4a) Reciprocity |S_HV / S_VH|", "dB", [("reciprocity_db", None)]),
        ("(4b) Reciprocity phase ∠(S_HV S_VH*)", "deg", [("reciprocity_phase_deg", None)]),
        ("(4c) Normalised asymmetry |S_HV − S_VH| / ‖S‖", "dB", [("asymmetry_db", None)]),
    ]
    for panel_index, (ax, (title, unit, series)) in enumerate(zip(axes, panels)):
        panel_flags = flagged if panel_index < 2 else flagged_crosspol
        for colour, (key, label) in zip(SERIES_COLOURS, series):
            band_plot(ax, x, metrics[key], valid, centre_index, colour, label, panel_flags)
        style_axis(ax, title, unit)
        if len(series) > 1:
            ax.legend(fontsize=7.5, frameon=False)
    if gating:
        for colour, (floor, label) in zip(
            SERIES_COLOURS, zip(crosspol_floor_db, ("Tx H floor", "Tx V floor"))
        ):
            axes[2].plot(x, floor, color=colour, linewidth=1.2, linestyle=":", label=label)
        axes[2].legend(fontsize=7, frameon=False)
    for level in (-20, -25, -30):
        axes[2].axhline(level, color=MUTED, linewidth=0.8, linestyle="--")
        axes[2].annotate(
            f"{level} dB",
            (x[-1], level),
            fontsize=7,
            color=MUTED,
            xytext=(2, 2),
            textcoords="offset points",
            ha="right",
        )

    ax = axes[6]
    names = ["k₁ HH+VV (odd)", "k₂ HH−VV (even)", "k₃ HV+VH", "k₄ j(HV−VH) (antisym.)"]
    for component, (colour, name) in enumerate(zip(SERIES_COLOURS, names)):
        band_plot(
            ax,
            x,
            metrics["pauli_db"][..., component],
            valid,
            centre_index,
            colour,
            name,
            flagged_crosspol,
        )
    style_axis(ax, "(5) Bistatic Pauli power fractions", "fraction of span (dB)")
    ax.legend(fontsize=7, frameon=False)

    ax = axes[7]
    band_plot(
        ax,
        x,
        metrics["infidelity_sphere_db"],
        valid,
        centre_index,
        SERIES_COLOURS[0],
        None,
        flagged_crosspol,
    )
    style_axis(ax, "(6) Infidelity to the sphere, 1 − F", "dB")

    ax = axes[8]
    band_plot(
        ax,
        x,
        metrics["alpha_deg"],
        valid,
        centre_index,
        SERIES_COLOURS[0],
        "single-target α (deg)",
        flagged_crosspol,
    )
    style_axis(ax, "(7) Cloude α, and entropy over the band", "α (deg)")
    ax_entropy = ax.inset_axes([0.55, 0.55, 0.42, 0.38])
    ax_entropy.plot(
        x,
        metrics["entropy_bistatic"],
        color=SERIES_COLOURS[1],
        marker="o",
        markersize=3,
        linewidth=1.5,
        label="H, 4×4 (bistatic)",
    )
    ax_entropy.plot(
        x,
        metrics["entropy_monostatic"],
        color=SERIES_COLOURS[2],
        marker="o",
        markersize=3,
        linewidth=1.5,
        label="H, 3×3 (symmetrised)",
    )
    ax_entropy.tick_params(labelsize=6.5, colors=MUTED)
    ax_entropy.grid(True, color=GRID, linewidth=0.5)
    ax_entropy.legend(fontsize=6, frameon=False)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")

    for ax in axes:
        ax.set_xlabel("co-polar spacing d (λ₀)", fontsize=8.5, color=MUTED)
        add_bistatic_axis(ax, x, wavelength, target_distance)
    f_lo, f_hi = frequencies[valid][[0, -1]] / 1e9
    title = (
        f"Boresight sphere, {'gated' if gating else 'ungated'}. "
        f"Line: {frequencies[centre_index] / 1e9:.2f} GHz; "
        f"band: min–max over {f_lo:.2f}–{f_hi:.2f} GHz."
    )
    if gating:
        title += (
            f"\nHollow markers: co-pol gate contaminated (< {MIN_GATE_SNR_DB:.0f} dB above the "
            f"guard band, or peak shifted > {1e3 * MAX_PATH_DEVIATION_M:.0f} mm); in panels using "
            f"cross-pol channels also cross-pol < {MIN_CROSSPOL_MARGIN_DB:.0f} dB above its floor "
            f"(dotted in (3))"
        )
    fig.suptitle(
        title,
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(os.path.join(out_dir, "signatures_vs_spacing.png"), dpi=150)
    plt.close(fig)


def linear_term_study(out_dir, wavelength, spacing_lambda):
    """
    Geometry-only illustration of the linear (plane-wave) term and of what survives compensating
    it, for the cluster's cross-polar pairs and a point target off boresight in the y-z plane.
    """
    frequency = np.array([SPEED_OF_LIGHT / wavelength])
    angles = np.linspace(0, 30, 301)
    directions = np.stack(
        [np.zeros_like(angles), np.sin(np.radians(angles)), np.cos(np.radians(angles))], axis=-1
    )

    def channel_phase(positions, rx, tx, model, target_range=None):
        tx_position, rx_position = pair_positions(positions, rx, tx)
        if model == "plane":
            return np.array(
                [plane_wave_phase(frequency, tx_position, rx_position, r)[0] for r in directions]
            )
        points = target_range * directions
        return spherical_wave_phase(frequency, tx_position, rx_position, points)[..., 0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    rows = []
    for colour, spacing in zip(
        SERIES_COLOURS, [spacing_lambda[0], spacing_lambda[3], spacing_lambda[-1]]
    ):
        positions = cluster_positions(spacing * wavelength)
        hv = channel_phase(positions, 0, 1, "plane")
        vh = channel_phase(positions, 1, 0, "plane")
        delta = np.degrees(np.unwrap(np.angle(hv / vh)))
        axes[0].plot(angles, delta, color=colour, linewidth=2, label=f"d = {spacing:.2f}λ")
        # A 45-deg dihedral (S_HV = S_VH) seen through the uncompensated linear term
        dihedral = np.zeros((len(angles), 2, 2), dtype=complex)
        dihedral[:, 0, 1], dihedral[:, 1, 0] = hv, vh
        axes[1].plot(
            angles,
            db10(pauli_power_fractions(dihedral)[:, 3]),
            color=colour,
            linewidth=2,
            label=f"d = {spacing:.2f}λ",
        )
        for angle in (1, 2, 5, 10, 20):
            index = np.argmin(np.abs(angles - angle))
            rows.append(
                (spacing, angle, delta[index], 100 * pauli_power_fractions(dihedral)[index, 3])
            )
    style_axis(axes[0], "Linear term: ∠S_HV − ∠S_VH, uncompensated", "deg")
    style_axis(axes[1], "45° dihedral: antisymmetric power, uncompensated", "k₄ fraction (dB)")
    axes[1].set_ylim(-40, 1)

    positions = cluster_positions(spacing_lambda[-1] * wavelength)
    fresnel_rows = []
    for colour, target_range in zip(SERIES_COLOURS, (0.5, 1.0, 3.0)):
        exact = {}
        for rx, tx in ((0, 0), (1, 1), (0, 1), (1, 0)):
            # Plane-wave compensation applied to the exact spherical-wave channel phase
            exact[(rx, tx)] = channel_phase(
                positions, rx, tx, "spherical", target_range
            ) / channel_phase(positions, rx, tx, "plane")
        residual = np.degrees(np.angle(exact[(0, 1)] / exact[(1, 0)]))
        residual_co = np.degrees(np.angle(exact[(0, 1)] / exact[(0, 0)]))
        axes[2].plot(
            angles, residual, color=colour, linewidth=2, label=f"HV−VH, R = {target_range:g} m"
        )
        axes[2].plot(
            angles,
            residual_co,
            color=colour,
            linewidth=1.2,
            linestyle="--",
            label=f"HV−HH, R = {target_range:g} m",
        )
        for angle in (5, 10, 20, 30):
            index = np.argmin(np.abs(angles - angle))
            fresnel_rows.append((target_range, angle, residual[index], residual_co[index]))
    style_axis(
        axes[2], f"Residual after plane-wave compensation, d = {spacing_lambda[-1]:.2f}λ", "deg"
    )
    for ax in axes:
        ax.set_xlabel("target elevation off boresight θ_y (deg)", fontsize=8.5, color=MUTED)
        ax.legend(fontsize=7, frameon=False)
    fig.savefig(os.path.join(out_dir, "linear_term_study.png"), dpi=150)
    plt.close(fig)
    return rows, fresnel_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("sweep", help=".npz written by pomaar.simulator.sbr_results")
    parser.add_argument("--out", default="data/processed/bistaticity")
    parser.add_argument("--wavelength", type=float, default=3.795e-3, help="λ₀ (m)")
    parser.add_argument("--target-range", type=float, default=3.0, help="Sphere centre (m)")
    parser.add_argument("--target-radius", type=float, default=0.1, help="Sphere radius (m)")
    parser.add_argument(
        "--feed-path",
        type=float,
        default=0.0505,
        help="Two-way feed path (m): 2 x 20 mm WR-12 waveguide at 79 GHz",
    )
    parser.add_argument(
        "--gate",
        choices=("auto", "on", "off"),
        default="auto",
        help="Range gating; auto gates only total-field data without a background",
    )
    parser.add_argument("--gate-half-width", type=float, default=0.06, help="Gate +- (m)")
    parser.add_argument(
        "--background",
        default=None,
        help="Same sweep without the target (.npz), subtracted coherently",
    )
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    processor = MimoPolarimetryProcessor.from_npz(args.sweep, args.background)
    raw = processor.scattering_matrices.copy()
    frequencies = processor.frequencies_hz
    spacing_lambda = processor.sweep_values
    setup_sweep = processor.metadata.get("setup_sweep", "")
    target_only = setup_sweep.endswith("_Scattered") or args.background is not None
    gating = args.gate == "on" or (args.gate == "auto" and not target_only)
    print(
        f"{processor.metadata.get('design')} ({setup_sweep}): {raw.shape[0]} spacings, "
        f"{len(frequencies)} points, {frequencies[0] / 1e9:.2f}-{frequencies[-1] / 1e9:.2f} GHz; "
        f"gating {'on' if gating else 'off'}"
        + (f"; background {args.background} subtracted" if args.background else "")
    )

    # Step 1: locate the sphere's specular return (common to all four channels)
    target_front = np.array([0.0, 0.0, args.target_range - args.target_radius])
    positions = [cluster_positions(spacing * args.wavelength) for spacing in spacing_lambda]
    target_paths = np.zeros(len(spacing_lambda))
    geometric_paths = np.zeros((len(spacing_lambda), 2, 2))
    print("spacing  expected path  found path  (m)")
    for index, spacing in enumerate(spacing_lambda):
        for rx in (0, 1):
            for tx in (0, 1):
                geometric_paths[index, rx, tx] = bistatic_path(
                    *pair_positions(positions[index], rx, tx), target_front
                )
        expected = geometric_paths[index].mean() + args.feed_path
        target_paths[index], _ = locate_target_path(raw[index], frequencies, expected, 0.03)
        print(f"{spacing:6.2f}λ  {expected:12.4f}  {target_paths[index]:10.4f}")

    # Step 1b: gate contamination: target peak inside the gate versus the guard band beside it.
    # Always reported; it only flags points when gating is on (otherwise there is nothing to leak).
    gate_peak = np.zeros((len(spacing_lambda), 2, 2))
    guard_peak = np.zeros((len(spacing_lambda), 2, 2))
    for index in range(len(spacing_lambda)):
        channels = np.moveaxis(raw[index], 0, -1)
        gate_peak[index] = delay_window_peak(
            channels,
            frequencies,
            target_paths[index],
            (-args.gate_half_width, args.gate_half_width),
            window=GATE_WINDOW,
        )
        guard_peak[index] = delay_window_peak(
            channels, frequencies, target_paths[index], GUARD_OFFSETS, window=GATE_WINDOW
        )
    gate_snr_db = db20(gate_peak / guard_peak)
    path_deviation = target_paths - np.median(target_paths)
    # Cross-pol leakage below this cannot be resolved: guard level of the cross-pol channel
    # relative to the co-pol target peak of the same transmit polarisation
    crosspol_floor_db = (
        db20(guard_peak[:, 1, 0] / gate_peak[:, 0, 0]),
        db20(guard_peak[:, 0, 1] / gate_peak[:, 1, 1]),
    )
    crosspol_margin_db = (
        db20(gate_peak[:, 1, 0] / gate_peak[:, 0, 0]) - crosspol_floor_db[0],
        db20(gate_peak[:, 0, 1] / gate_peak[:, 1, 1]) - crosspol_floor_db[1],
    )
    flagged = np.zeros(len(spacing_lambda), dtype=bool)
    flagged_crosspol = np.zeros(len(spacing_lambda), dtype=bool)
    if gating:
        flagged = (
            (gate_snr_db[:, 0, 0] < MIN_GATE_SNR_DB)
            | (gate_snr_db[:, 1, 1] < MIN_GATE_SNR_DB)
            | (np.abs(path_deviation) > MAX_PATH_DEVIATION_M)
        )
        flagged_crosspol = flagged | (np.minimum(*crosspol_margin_db) < MIN_CROSSPOL_MARGIN_DB)
    print("gate SNR (dB) HH/VV/HV/VH per spacing:")
    for index, spacing in enumerate(spacing_lambda):
        snr = gate_snr_db[index]
        print(
            f"  {spacing:5.2f}λ  {snr[0, 0]:6.1f} {snr[1, 1]:6.1f} {snr[0, 1]:6.1f} "
            f"{snr[1, 0]:6.1f}  peak shift {1e3 * path_deviation[index]:+5.1f} mm  x-pol margin "
            f"{crosspol_margin_db[0][index]:5.1f}/{crosspol_margin_db[1][index]:5.1f} dB"
            f"{'  <- co-pol contaminated' if flagged[index] else ''}"
            f"{'  <- x-pol unresolved' if flagged_crosspol[index] > flagged[index] else ''}"
        )
    if gating:
        processor.gate(target_paths, args.gate_half_width, window=GATE_WINDOW)
    valid = processor.valid_band
    centre_index = np.argmin(np.abs(frequencies - SPEED_OF_LIGHT / args.wavelength))

    # Step 2: geometric compensation of each channel (spherical wave, focus on the sphere front)
    uncompensated = copy.deepcopy(processor)
    processor.compensate_geometry(positions, target_front)
    channel_path_spread = geometric_paths.max(axis=(1, 2)) - geometric_paths.min(axis=(1, 2))
    print(
        "max inter-channel geometric path difference (um):", np.round(channel_path_spread * 1e6, 3)
    )

    # Step 3: signatures (on the compensated matrices; for this square cluster the ratios are
    # identical to the uncompensated ones)
    metrics = plot_metrics(processor.signatures(REFERENCE_SPHERE))
    metrics_uncompensated = plot_metrics(uncompensated.signatures(REFERENCE_SPHERE))
    change = np.max(np.abs(metrics["copol_phase_deg"] - metrics_uncompensated["copol_phase_deg"]))
    print(f"compensation changes the co-pol phase difference by at most {change:.2e} deg")

    processed = processor.scattering_matrices
    target_distance = args.target_range - args.target_radius
    plot_validation(
        args.out,
        raw,
        processed,
        frequencies,
        spacing_lambda,
        target_paths,
        args.gate_half_width,
        valid,
        centre_index,
        gate_snr_db,
        flagged,
        gating,
    )
    plot_signatures(
        args.out,
        metrics,
        spacing_lambda,
        valid,
        centre_index,
        args.wavelength,
        target_distance,
        frequencies,
        flagged,
        flagged_crosspol,
        crosspol_floor_db,
        gating,
    )
    linear_rows, fresnel_rows = linear_term_study(args.out, args.wavelength, spacing_lambda)

    with open(os.path.join(args.out, "signatures_vs_spacing.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        scalar_keys = [k for k, v in metrics.items() if v.ndim == 2]
        writer.writerow(
            [
                "spacing_lambda",
                "beta_deg",
                "abs_s_hh_db",
                "gated",
                "contaminated",
                "crosspol_unresolved",
                "gate_snr_hh_db",
                "gate_snr_vv_db",
                "xpol_floor_txh_db",
                "xpol_floor_txv_db",
            ]
            + [f"{k}@f0" for k in scalar_keys]
            + [f"{k}_band_{e}" for k in scalar_keys for e in ("min", "max")]
            + [
                "pauli_k1_db@f0",
                "pauli_k2_db@f0",
                "pauli_k3_db@f0",
                "pauli_k4_db@f0",
                "entropy_bistatic",
                "entropy_monostatic",
            ]
        )
        for index, spacing in enumerate(spacing_lambda):
            beta = np.degrees(2 * np.arctan(spacing * args.wavelength / (2 * target_distance)))
            writer.writerow(
                [
                    spacing,
                    beta,
                    db20(processed[index, centre_index, 0, 0]),
                    gating,
                    bool(flagged[index]),
                    bool(flagged_crosspol[index]),
                    gate_snr_db[index, 0, 0],
                    gate_snr_db[index, 1, 1],
                    crosspol_floor_db[0][index],
                    crosspol_floor_db[1][index],
                ]
                + [metrics[k][index, centre_index] for k in scalar_keys]
                + [f(metrics[k][index, valid]) for k in scalar_keys for f in (np.min, np.max)]
                + list(metrics["pauli_db"][index, centre_index])
                + [metrics["entropy_bistatic"][index], metrics["entropy_monostatic"][index]]
            )
    with open(os.path.join(args.out, "linear_term.csv"), "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["spacing_lambda", "theta_y_deg", "dphi_hv_vh_deg", "dihedral45_antisym_percent"]
        )
        writer.writerows(linear_rows)
        writer.writerow([])
        writer.writerow(
            ["target_range_m", "theta_y_deg", "residual_hv_vh_deg", "residual_hv_hh_deg"]
        )
        writer.writerows(fresnel_rows)
    print(f"Figures and CSVs written to {args.out}")


if __name__ == "__main__":
    main()
