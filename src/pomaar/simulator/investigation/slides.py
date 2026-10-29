#!/usr/bin/env python3
"""
Presentation figures of the parametric-beam bistaticity investigation (host side, NumPy only):
one message per figure, slide-sized, from the exports in <exports> and the horn-cluster export in
data/interim/bistaticity/:

    python src/pomaar/simulator/investigation/slides.py [--out <folder>]

Writes to reports/bistaticity-analysis/parametric_beams/ by default:
* `01_scene.png`: the multi-receiver ladder and the split pairs (schematic, not to scale);
* `02_horn_vs_ideal.png`: LDR and cross-co phase of a boresight sphere, horn cluster with
  embedded patterns against parametric beams;
* `03_bistatic_patterns.png`: the span against beta per target (the targets' bistatic patterns);
* `04_lobe_scaling.png`: the trihedral's span over size and range, against beta and beta / beta_0;
* `05_copol_lobe_limited.png`: span, Z_DR and co-pol phase change against beta / beta_0;
* `06_symmetry_order.png`: emerging channels and co-pol phase change on log-log axes;
* `07_cross_pol_genuine_vs_floor.png`: LDR_H and the cross-co phase, genuine against floor;
* `08_cloude_pottier.png`: LDR_H against entropy and mean alpha of the same targets;
* `09_cross_pol_size_range.png`: the trihedral's cross-pol over size and range;
* `10_near_field.png`: a dihedral inside its far-field distance against one outside it;
* `11_split_pairs.png`: Z_DR from split pairs minus a co-located pair;
* `12_rx_aim.png`: a copy of <results>/rx_aim/rx_aim_scalar_fit.png.

Every curve is the centre frequency, every band the min-max over 74-84 GHz; the receivers are
aimed at the target (rxAim 1). Like the other results, the figures show the data only: no
reference lines, thresholds or expected values (the slope and departure numbers are in the log).
"""

import argparse
import json
from pathlib import Path
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.config import INTERIM_DATA_DIR, REPORTS_DIR  # noqa: E402
from pomaar.simulator.bistaticity_analysis import (  # noqa: E402
    LABELS,
    band_limits,
    pool_receivers,
    wrap_deg,
)
from pomaar.simulator.investigation.pairs import load as load_pairs  # noqa: E402
from pomaar.simulator.investigation.project import EXPORTS_DIR, RESULTS_DIR  # noqa: E402
from pomaar.simulator.investigation.scaling import first_null  # noqa: E402
from pomaar.simulator.polarimetry_processor import CanonicalTargetProcessor  # noqa: E402

OUT_DIR = REPORTS_DIR / "bistaticity-analysis" / "parametric_beams"
HORN_SPHERE = INTERIM_DATA_DIR / "bistaticity" / "horn_cluster_sphere_boresight_3m.npz"
INK, MUTED, GRID = "#222222", "#52514e", "#e4e4e4"
WIDE, TALL = 13.33, 5.4
# one colour per target throughout (the colour follows the target, never its rank)
TARGETS = {
    "MultiRxTrihedral_rxAim1_targetRoll0": ("Trihedral, roll 0°", "#2a78d6"),
    "MultiRxTrihedral_rxAim1_targetRoll90": ("Trihedral, roll 90°", "#eb6834"),
    "MultiRxTrihedral_rxAim1_targetRoll45": ("Trihedral, roll 45°", "#1baf7a"),
    "MultiRxDihedral_rxAim1_targetRoll0": ("Dihedral, roll 0°", "#eda100"),
    "MultiRxDihedral_rxAim1_targetRoll90": ("Dihedral, roll 90°", "#e87ba4"),
    "MultiRxPlate_rxAim1": ("Plate", "#008300"),
    "MultiRxSphere_rxAim1": ("Sphere", "#6250d6"),
    "MultiRxDihedral_rxAim1_targetRoll45": ("Dihedral, roll 45°", "#e34948"),
    "MultiRxDihedral_rxAim1_targetRoll22.5": ("Dihedral, roll 22.5°", "#7a7a7a"),
}
# trihedral size and range: colour by size, line style by range
SCALING = {
    "MultiRxTrihedralR1_rxAim1": ("40 mm, 1 m", "#2a78d6", "--"),
    "MultiRxTrihedral_rxAim1_targetRoll0": ("40 mm, 4 m", "#2a78d6", "-"),
    "MultiRxTrihedralS80_rxAim1": ("80 mm, 4 m", "#1baf7a", "-"),
    "MultiRxTrihedralS80R16_rxAim1": ("80 mm, 16 m", "#1baf7a", ":"),
    "MultiRxTrihedralS160_rxAim1": ("160 mm, 4 m", "#6250d6", "-"),
}
NEAR_FIELD = {
    "MultiRxDihedral_rxAim1_targetRoll90": ("40 mm, far field from 1.7 m", "#e87ba4"),
    "MultiRxDihedralS160_rxAim1": ("160 mm, far field from 27 m", "#7a7a7a"),
}
PAIRS = {
    "PairsSphere": ("Sphere", "#6250d6"),
    "PairsTrihedral_targetRoll0": ("Trihedral, roll 0°", "#2a78d6"),
    "PairsTrihedral_targetRoll90": ("Trihedral, roll 90°", "#eb6834"),
    "PairsDihedral_targetRoll0": ("Dihedral, roll 0°", "#eda100"),
    "PairsDihedral_targetRoll22.5": ("Dihedral, roll 22.5°", "#7a7a7a"),
}
# a channel counts as emerged where it stands this far above its monostatic value
EMERGED_DB = 10.0
BETA = r"bistatic angle $\beta$ (deg)"
BETA_RATIO = r"$\beta / \beta_0$ ($\beta_0$: first null of the span)"

plt.rcParams.update(
    {
        "font.size": 13,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "legend.fontsize": 11.5,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.edgecolor": GRID,
        "axes.labelcolor": MUTED,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 200,
        "savefig.facecolor": "white",
    }
)

_cache = {}


def pooled(stem):
    """(features, angles, centre index) of an export, pooled by bistatic angle."""
    if stem not in _cache:
        features, angles, _, _ = pool_receivers(EXPORTS_DIR / f"{stem}.npz")
        centre = features["span_db"].shape[1] // 2
        _cache[stem] = features, angles, centre
    return _cache[stem]


def change(values, circular):
    """Each row minus the monostatic one (first angle), wrapped for phases."""
    difference = values - values[:1]
    return wrap_deg(difference) if circular else difference


def curve(ax, x, values, centre, colour, label, circular=False, style="-", band=True, mask=None):
    """Centre-frequency line with the band's min-max (per frequency values) where `mask`."""
    keep = np.ones(len(x), bool) if mask is None else mask
    if values.ndim == 2 and band:
        low, high = band_limits(values, centre, circular)
        ax.fill_between(x, np.where(keep, low, np.nan), np.where(keep, high, np.nan),
                        color=colour, alpha=0.13, lw=0)  # fmt: skip
    line = values[:, centre] if values.ndim == 2 else values
    ax.plot(x, np.where(keep, line, np.nan), color=colour, lw=2.2, ls=style, label=label,
            marker="o", markersize=3.5)  # fmt: skip


def finish(ax, title, ylabel, xlabel=BETA):
    ax.set_title(title, loc="left", color=INK)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)


def phase_axis(ax):
    ax.set_ylim(-180, 180)
    ax.set_yticks(np.arange(-180, 181, 90))


def save(fig, out, name):
    path = out / name
    fig.savefig(path)
    plt.close(fig)
    print(f"Wrote {path}")


def scene(out):
    """Schematic top views (x-y, not to scale): multi-receiver ladder and split pairs."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(WIDE, TALL))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.08, wspace=0.08)
    for ax in (ax1, ax2):
        ax.set_aspect("equal")
        ax.axis("off")
    target = np.array([10.0, 0.0])
    beta = np.radians(9.0)
    ys = [0.0, target[0] * np.tan(beta), target[0] * np.tan(2 * beta)]
    ax1.plot(*target, marker="D", ms=13, color=INK)
    ax1.annotate("target at range R\n(boresight, +x)", target, xytext=(0, -38),
                 textcoords="offset points", ha="center", color=INK)  # fmt: skip
    colours = ["#2a78d6", "#1baf7a", "#6250d6"]
    for k, (y, colour) in enumerate(zip(ys, colours)):
        ax1.plot([target[0], 0], [0, y], color=colour, lw=1.6, ls="--" if k else "-")
        ax1.plot(0, y, marker="s", ms=13, color=colour, mfc="white" if k else colour, mew=2.2)
        name = "Tx and Rx$_0$ (H+V)" if k == 0 else f"Rx$_{k}$ (H+V)"
        ax1.annotate(name, (0, y), xytext=(-12, 0), textcoords="offset points", ha="right",
                     va="center", color=INK)  # fmt: skip
        if k:
            arc = np.linspace(np.pi, np.pi - k * beta, 30)
            radius = 3.2 + 0.9 * k
            ax1.plot(target[0] + radius * np.cos(arc), radius * np.sin(arc), color=colour, lw=1.4)
            mid = arc[len(arc) // 2]
            ax1.annotate(r"$\beta$" if k == 1 else r"$2\beta$",
                         (target[0] + (radius + 0.5) * np.cos(mid), (radius + 0.5) * np.sin(mid)),
                         ha="center", va="center", color=colour)  # fmt: skip
    ax1.annotate("", (0, ys[-1] + 1.2), (0, -1.2),
                 arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.2))  # fmt: skip
    ax1.text(0.35, ys[-1] + 1.0, "y (H)", color=MUTED)
    ax1.set_xlim(-6.5, 12.5)
    ax1.set_ylim(-3.2, ys[-1] + 2.0)
    fig.text(0.01, 0.93, "(a) Multi-receiver ladder: $y_k = R\\,\\tan(k\\beta)$, one solve "
             "gives 0, $\\beta$, 2$\\beta$", color=INK, fontsize=14)  # fmt: skip

    # split pairs: TxV -b- RxV -g- TxH -b- RxH along y, target on boresight
    b, g = 1.6, 3.2
    positions = {"TxV": 0.0, "RxV": b, "TxH": b + g, "RxH": 2 * b + g}
    centre = (2 * b + g) / 2
    positions = {k: v - centre for k, v in positions.items()}
    target = np.array([10.0, 0.0])
    ax2.plot(*target, marker="D", ms=13, color=INK)
    for name, y in positions.items():
        colour = "#eb6834" if name.endswith("V") else "#2a78d6"
        filled = name.startswith("Tx")
        ax2.plot([target[0], 0], [0, y], color=colour, lw=1.4, ls="-" if filled else "--")
        ax2.plot(0, y, marker="s", ms=13, color=colour, mfc=colour if filled else "white", mew=2.2)
        ax2.annotate(name, (0, y), xytext=(-12, 0), textcoords="offset points", ha="right",
                     va="center", color=INK)  # fmt: skip
    for (start, stop), text in (
        (("TxV", "RxV"), "b"), (("RxV", "TxH"), "g"), (("TxH", "RxH"), "b")
    ):  # fmt: skip
        y0, y1 = positions[start], positions[stop]
        ax2.annotate("", (0.9, y1), (0.9, y0),
                     arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.1))  # fmt: skip
        ax2.text(1.2, (y0 + y1) / 2, text, va="center", color=MUTED)
    ax2.text(5.0, -4.2, "VV from TxV–RxV, HH from TxH–RxH: same pair angle $\\beta$,\n"
             "bisectors $(b+g)/R$ apart", ha="center", va="top", color=INK)  # fmt: skip
    ax2.set_xlim(-6.5, 12.5)
    ax2.set_ylim(-6.2, 5.5)
    fig.text(0.52, 0.93, "(b) Split $Z_\\mathrm{DR}$ pairs against a co-located pair", color=INK,
             fontsize=14)  # fmt: skip
    fig.text(0.01, 0.01, "Top views, not to scale. Parametric beams (60° × 60°, identical H and "
             "V, no pattern cross-pol), H along y and V along z.", color=MUTED,
             fontsize=11.5)  # fmt: skip
    save(fig, out, "01_scene.png")


def horn_sphere():
    """Horn cluster, sphere at 3 m: (features, co-pol bistatic angle per spacing, centre)."""
    processor = CanonicalTargetProcessor.from_npz(HORN_SPHERE)
    features = processor.signatures(matrices=processor.target_matrices(0.0, 0.0))
    metadata = json.loads(str(np.load(HORN_SPHERE)["metadata"]))
    nominal = metadata["nominal_variables"]
    wavelength = float(nominal["centreWavelength"].removesuffix("mm")) * 1e-3
    distance = float(nominal["targetRange"].removesuffix("meter"))
    spacing = np.asarray(processor.sweep_values, float) * wavelength
    beta = np.degrees(2 * np.arctan(spacing / (2 * distance)))
    return features, beta, features["span_db"].shape[1] // 2


def horn_vs_ideal(out):
    horn, horn_beta, horn_centre = horn_sphere()
    ideal, angles, centre = pooled("MultiRxSphere_rxAim1")
    show = angles <= 1.0
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(WIDE, TALL), constrained_layout=True)
    for features, x, c, colour, label, keep in (
        (horn, horn_beta, horn_centre, INK, "Horn cluster, embedded patterns (r = 100 mm, 3 m)",
         None),
        (ideal, angles, centre, "#6250d6", "Parametric beams (r = 20 mm, 4 m)", show),
    ):  # fmt: skip
        x_shown = x if keep is None else x[keep]
        for ax, key, circular in (
            (ax1, "ldr_txh_db", False),
            (ax2, "cross_co_phase_txh_deg", True),
        ):
            values = features[key] if keep is None else features[key][keep]
            curve(ax, x_shown, values, c, colour, label, circular)
    finish(ax1, LABELS["ldr_txh_db"], "dB")
    finish(ax2, LABELS["cross_co_phase_txh_deg"], "deg")
    phase_axis(ax2)
    fig.legend(*ax1.get_legend_handles_labels(), loc="outside lower center", ncol=2,
               frameon=False)  # fmt: skip
    fig.suptitle("Sphere on boresight, the same processing: horn cluster (bistatic angle of the "
                 "co-pol pairs) against ideal parametric beams", x=0.01, ha="left", color=INK,
                 fontsize=14)  # fmt: skip
    save(fig, out, "02_horn_vs_ideal.png")


def bistatic_patterns(out):
    fig, ax = plt.subplots(figsize=(WIDE, TALL), constrained_layout=True)
    for stem in (
        "MultiRxSphere_rxAim1",
        "MultiRxPlate_rxAim1",
        "MultiRxDihedral_rxAim1_targetRoll0",
        "MultiRxDihedral_rxAim1_targetRoll90",
        "MultiRxTrihedral_rxAim1_targetRoll0",
    ):
        features, angles, centre = pooled(stem)
        label, colour = TARGETS[stem]
        curve(ax, angles, change(features["span_db"], False), centre, colour, label)
    finish(ax, f"Span {LABELS['span_db']} relative to monostatic (all targets 40 mm across, "
           "at 4 m)", "dB")  # fmt: skip
    ax.set_ylim(-45, 3)
    ax.legend(frameon=False, loc="lower left")
    save(fig, out, "03_bistatic_patterns.png")


def lobe_scaling(out):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(WIDE, TALL), constrained_layout=True)
    for stem, (label, colour, style) in SCALING.items():
        features, angles, centre = pooled(stem)
        span = change(features["span_db"], False)
        null = first_null(angles, span[:, centre])
        curve(ax1, angles, span, centre, colour, label, style=style, band=False)
        curve(ax2, angles / null, span, centre, colour,
              rf"{label}: $\beta_0$ = {null:.2f}°", style=style, band=False)  # fmt: skip
    finish(ax1, "Trihedral (roll 0°): span relative to monostatic", "dB")
    finish(ax2, "The same against the first null", "dB", BETA_RATIO)
    for ax in (ax1, ax2):
        ax.set_ylim(-45, 3)
    ax1.set_xlim(0, 12)
    ax2.set_xlim(0, 2)
    ax2.legend(frameon=False, loc="lower left")
    save(fig, out, "04_lobe_scaling.png")


def copol_lobe_limited(out):
    stems = [
        "MultiRxTrihedral_rxAim1_targetRoll0",
        "MultiRxTrihedral_rxAim1_targetRoll90",
        "MultiRxTrihedral_rxAim1_targetRoll45",
        "MultiRxDihedral_rxAim1_targetRoll0",
        "MultiRxDihedral_rxAim1_targetRoll90",
        "MultiRxPlate_rxAim1",
    ]
    fig, axes = plt.subplots(1, 3, figsize=(WIDE, TALL), constrained_layout=True)
    for stem in stems:
        features, angles, centre = pooled(stem)
        label, colour = TARGETS[stem]
        null = first_null(angles, features["span_db"][:, centre])
        x = angles / null
        keep = x <= 1.3
        for ax, key, circular in zip(
            axes, ("span_db", "zdr_db", "hh_vv_phase_deg"), (False, False, True)
        ):
            values = change(features[key], circular)[keep]
            curve(ax, x[keep], values, centre, colour, label, circular, band=False)
    finish(axes[0], f"Span {LABELS['span_db']}", "dB, change from monostatic", BETA_RATIO)
    finish(axes[1], LABELS["zdr_db"], "dB, change from monostatic", BETA_RATIO)
    finish(axes[2], LABELS["hh_vv_phase_deg"], "deg, change from monostatic", BETA_RATIO)
    axes[0].set_ylim(-40, 3)
    axes[1].set_ylim(-6, 6)
    axes[2].set_ylim(-60, 60)
    axes[0].legend(frameon=False, loc="lower left", fontsize=10.5)
    save(fig, out, "05_copol_lobe_limited.png")


def symmetry_order(out):
    """Log-log: emerging channels (power re monostatic span) and co-pol phase change."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(WIDE, TALL), constrained_layout=True)
    emerging = [
        ("MultiRxTrihedral_rxAim1_targetRoll0", "power_vh_db", r"$|S_\mathrm{VH}|^2$"),
        ("MultiRxTrihedral_rxAim1_targetRoll45", "power_vh_db", r"$|S_\mathrm{VH}|^2$"),
        ("MultiRxDihedral_rxAim1_targetRoll45", ("pauli_fractions_db", 1), r"$|k_2|^2$"),
        ("MultiRxDihedral_rxAim1_targetRoll22.5", ("pauli_fractions_db", 3), r"$|k_4|^2$"),
    ]
    for stem, key, name in emerging:
        features, angles, centre = pooled(stem)
        label, colour = TARGETS[stem]
        span0 = features["span_db"][0, centre]
        if isinstance(key, tuple):
            values = features[key[0]][:, centre, key[1]] + features["span_db"][:, centre]
        else:
            values = features[key][:, centre]
        values = values - span0
        keep = (angles > 0) & (angles <= 4.0) & (values > values[0] + EMERGED_DB)
        ax1.plot(angles[keep], values[keep], color=colour, lw=2.2, marker="o", markersize=4,
                 label=f"{label}: {name}")  # fmt: skip
    ax1.set_xscale("log")
    finish(ax1, "Channels that are zero monostatically", "dB re monostatic span")
    ax1.legend(frameon=False, loc="lower right")
    for stem in (
        "MultiRxTrihedral_rxAim1_targetRoll0",
        "MultiRxTrihedral_rxAim1_targetRoll90",
        "MultiRxTrihedral_rxAim1_targetRoll45",
    ):
        features, angles, centre = pooled(stem)
        label, colour = TARGETS[stem]
        values = np.abs(change(features["hh_vv_phase_deg"], True)[:, centre])
        keep = (angles > 0) & (angles <= 4.0)
        ax2.plot(angles[keep], values[keep], color=colour, lw=2.2, marker="o", markersize=4,
                 label=label)  # fmt: skip
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    finish(ax2, f"Co-pol phase {LABELS['hh_vv_phase_deg']}: |change from monostatic|", "deg")
    ax2.legend(frameon=False, loc="lower right")
    for ax in (ax1, ax2):
        ax.set_xticks([0.25, 0.5, 1, 2, 4], ["0.25", "0.5", "1", "2", "4"])
        ax.minorticks_off()
    ax2.set_yticks([0.01, 0.1, 1, 10], ["0.01", "0.1", "1", "10"])
    save(fig, out, "06_symmetry_order.png")


def broken(phase_deg):
    """A wrapped phase with NaN at its wraps, so that a line does not jump across the axis."""
    phase = np.asarray(phase_deg, float).copy()
    jumps = np.flatnonzero(np.abs(np.diff(phase)) > 180.0)
    return np.insert(phase, jumps + 1, np.nan), jumps


def cross_pol_genuine_vs_floor(out):
    """
    The trihedral (roll 0) at bistatic angles 0, 1 and 2 deg (trust check T4) at 4, 8 and 16
    rays per wavelength: LDR_H and the cross-co phase against frequency per receiver.
    """
    runs = {4: "t1_trihedral_1deg", 8: "t4_trihedral_1deg_rd8", 16: "t4_trihedral_1deg_rd16"}
    receivers = {
        "Rx_0/Tx": ("0° (monostatic)", "#2a78d6"),
        "Rx_1/Tx": ("1°", "#1baf7a"),
        "Rx_2/Tx": ("2°", "#6250d6"),
    }
    styles = {4: "-", 8: "--", 16: ":"}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(WIDE, TALL), constrained_layout=True)
    for density, stem in runs.items():
        path = EXPORTS_DIR / f"{stem}.npz"
        frequencies = np.load(path)["frequencies_hz"] / 1e9
        for quad, (name, colour) in receivers.items():
            processor = CanonicalTargetProcessor.from_npz(path, quad=quad)
            features = processor.signatures(matrices=processor.target_matrices(0.0, 0.0))
            label = rf"$\beta$ = {name}, {density} rays/$\lambda$"
            style = dict(color=colour, ls=styles[density], lw=2.0)
            ax1.plot(frequencies, features["ldr_txh_db"][0], label=label, **style)
            phase, jumps = broken(features["cross_co_phase_txh_deg"][0])
            ax2.plot(np.insert(frequencies, jumps + 1, np.nan), phase, **style)
    finish(ax1, f"Trihedral (roll 0°): {LABELS['ldr_txh_db']}", "dB", "frequency (GHz)")
    finish(ax2, LABELS["cross_co_phase_txh_deg"], "deg", "frequency (GHz)")
    phase_axis(ax2)
    fig.legend(*ax1.get_legend_handles_labels(), loc="outside right center", frameon=False,
               fontsize=10.5)  # fmt: skip
    save(fig, out, "07_cross_pol_genuine_vs_floor.png")


def cloude_pottier(out):
    stems = ["MultiRxTrihedral_rxAim1_targetRoll0", "MultiRxTrihedral_rxAim1_targetRoll45"]
    fig, axes = plt.subplots(1, 3, figsize=(WIDE, TALL), constrained_layout=True)
    for stem in stems:
        features, angles, centre = pooled(stem)
        label, colour = TARGETS[stem]
        keep = angles <= 10.0
        curve(axes[0], angles[keep], features["ldr_txh_db"][keep], centre, colour, label)
        curve(axes[1], angles[keep], features["entropy"][keep], centre, colour, label)
        curve(axes[2], angles[keep], features["mean_alpha_deg"][keep], centre, colour, label)
    finish(axes[0], LABELS["ldr_txh_db"], "dB")
    finish(axes[1], "Entropy $H$ (band-averaged)", "")
    finish(axes[2], r"Mean $\bar{\alpha}$ (band-averaged)", "deg")
    axes[1].set_ylim(0, 1)
    axes[2].set_ylim(0, 90)
    axes[0].legend(frameon=False, loc="lower right")
    save(fig, out, "08_cloude_pottier.png")


def cross_pol_size_range(out):
    fig, ax = plt.subplots(figsize=(WIDE, TALL), constrained_layout=True)
    for stem, (label, colour, style) in SCALING.items():
        features, angles, centre = pooled(stem)
        values = features["power_hv_db"][:, centre] - features["span_db"][0, centre]
        keep = (angles > 0) & (values > values[0] + EMERGED_DB)
        ax.plot(angles[keep], values[keep], color=colour, ls=style, lw=2.2, marker="o",
                markersize=4, label=label)  # fmt: skip
    ax.set_xscale("log")
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8], ["0.25", "0.5", "1", "2", "4", "8"])
    ax.minorticks_off()
    finish(ax, f"Trihedral (roll 0°): {LABELS['power_hv_db']} relative to the monostatic span",
           "dB")  # fmt: skip
    ax.legend(frameon=False, loc="lower right")
    save(fig, out, "09_cross_pol_size_range.png")


def near_field(out):
    fig, axes = plt.subplots(1, 3, figsize=(WIDE, TALL), constrained_layout=True)
    for stem, (label, colour) in NEAR_FIELD.items():
        features, angles, centre = pooled(stem)
        keep = angles <= 6.0
        for ax, key, circular in zip(
            axes, ("span_db", "zdr_db", "hh_vv_phase_deg"), (False, False, True)
        ):
            values = change(features[key], circular)[keep]
            curve(ax, angles[keep], values, centre, colour, label, circular)
    finish(axes[0], f"Span {LABELS['span_db']}", "dB, change from monostatic")
    finish(axes[1], LABELS["zdr_db"], "dB, change from monostatic")
    finish(axes[2], LABELS["hh_vv_phase_deg"], "deg, change from monostatic")
    axes[1].set_ylim(-3, 3)
    axes[2].set_ylim(-30, 30)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=2,
               frameon=False, title="Dihedral (roll 90°) at 4 m", title_fontsize=11.5)  # fmt: skip
    save(fig, out, "10_near_field.png")


def split_pairs(out):
    fig, axes = plt.subplots(1, 2, figsize=(WIDE, TALL), constrained_layout=True, sharey=True)
    for design, (label, colour) in PAIRS.items():
        for ax, gap in zip(axes, ("0", "2")):
            stem = (
                design.replace("_", f"_gapRatio{gap}_", 1)
                if "_" in design
                else (f"{design}_gapRatio{gap}")
            )
            run = load_pairs(EXPORTS_DIR / f"{stem}.npz")
            centre = len(run["frequencies_hz"]) // 2
            curve(ax, run["pair_angle_deg"], run["zdr_split_minus_ref_db"], centre, colour,
                  label)  # fmt: skip
    finish(axes[0], r"$\Delta Z_\mathrm{DR}$, split minus co-located, gap $g = 0$",
           "dB", r"pair bistatic angle $\beta$ (deg)")  # fmt: skip
    finish(axes[1], r"The same with gap $g = 2b$", "", r"pair bistatic angle $\beta$ (deg)")
    axes[0].set_ylim(-10, 10)
    axes[0].legend(frameon=False, loc="lower left")
    save(fig, out, "11_split_pairs.png")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default=str(OUT_DIR))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    scene(out)
    horn_vs_ideal(out)
    bistatic_patterns(out)
    lobe_scaling(out)
    copol_lobe_limited(out)
    symmetry_order(out)
    cross_pol_genuine_vs_floor(out)
    cloude_pottier(out)
    cross_pol_size_range(out)
    near_field(out)
    split_pairs(out)
    shutil.copyfile(RESULTS_DIR / "rx_aim" / "rx_aim_scalar_fit.png", out / "12_rx_aim.png")
    print(f"Wrote {out / '12_rx_aim.png'} (copy)")


if __name__ == "__main__":
    main()
