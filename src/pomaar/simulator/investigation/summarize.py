#!/usr/bin/env python3
"""
Cross-target summary of the multi-receiver bistaticity sweep (host side, NumPy only): every
export in <sandbox>/exports matching the given stems is pooled by bistatic angle
(`bistaticity_analysis.pool_receivers`) and compared, at the centre frequency with the min-max
band over the sweep:

    python src/pomaar/simulator/investigation/summarize.py --aim 1 [--out <folder>]

Figures (in <sandbox>/summary/):
* `summary_signatures_aim<a>.png`: selected signatures against the bistatic angle, one series
  per design and roll;
* `summary_change_aim<a>.png`: the signatures that are well defined at the monostatic point
  (not set by the numerical cross-pol floor) minus their monostatic value, i.e. the change that
  the bistatic geometry alone causes, per target;
* `summary_aim_difference.png`: rxAim 1 minus rxAim 0 (receivers aimed at the target against
  receivers on boresight), i.e. the antenna-pattern and geometric-squint part;
* `summary_aim<a>.csv`: the centre-frequency values and band limits behind the figures.
"""

import argparse
import csv
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.simulator.bistaticity_analysis import (  # noqa: E402
    INK,
    MUTED,
    PAULI_NAMES,
    band_limits,
    pool_receivers,
    style_axis,
    wrap_deg,
)
from pomaar.simulator.investigation.sandbox import SANDBOX_DIR  # noqa: E402

EXPORTS = SANDBOX_DIR / "exports"
# (title, unit, feature key or (key, index), circular, domain)
PANELS = (
    [
        ("Span ‖S‖²", "dB", "span_db", False, None),
        ("Co-pol ratio |S_VV / S_HH|", "dB", "copol_ratio_db", False, None),
        ("Co-pol phase difference ∠(S_VV S_HH*)", "deg", "copol_phase_difference_deg", True, None),
        ("Co/cross ratio, Tx H: |S_HH / S_VH|", "dB", "co_minus_cross_txh_db", False, None),
        ("Co/cross phase, Tx H: ∠(S_HH S_VH*)", "deg", "co_cross_phase_txh_deg", True, None),
        ("Co/cross ratio, Tx V: |S_VV / S_HV|", "dB", "co_minus_cross_txv_db", False, None),
        ("Cross-pol ratio |S_HV / S_VH|", "dB", "reciprocity_db", False, None),
        ("Normalized asymmetry |S_HV − S_VH| / ‖S‖", "dB", "asymmetry_db", False, None),
    ]
    + [
        (f"Pauli fraction {name}", "dB", ("pauli_fractions_db", i), False, None)
        for i, name in enumerate(PAULI_NAMES)
    ]
    + [
        ("Cloude–Pottier entropy H (band-averaged)", "", "entropy", False, (0.0, 1.0)),
        ("Cloude–Pottier mean α (band-averaged)", "deg", "mean_alpha_deg", False, (0.0, 90.0)),
    ]
)
COLOURS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8f5bd6", "#d6455d", "#3fb6c9", "#7a7a7a",
]  # fmt: skip


def label_of(stem):
    """'MultiRxDihedral_rxAim1_targetRoll22.5' -> 'Dihedral, roll 22.5°'."""
    design = re.sub(r"^MultiRx", "", stem.split("_")[0])
    roll = re.search(r"targetRoll([\d.]+)", stem)
    return f"{design}, roll {float(roll[1]):g}°" if roll else design


def values_of(features, key):
    return features[key[0]][..., key[1]] if isinstance(key, tuple) else features[key]


def load(stems):
    pooled = {}
    for stem in stems:
        features, angles, receivers, _ = pool_receivers(EXPORTS / f"{stem}.npz")
        pooled[stem] = (features, angles)
    return pooled


def draw(ax, x, values, centre, colour, label, circular):
    if values.ndim == 1:
        ax.plot(x, values, color=colour, lw=1.6, marker="o", markersize=2.5, label=label)
        return
    low, high = band_limits(values, centre, circular)
    line = values[:, centre]
    for shift in (-360.0, 0.0, 360.0) if circular else (0.0,):
        ax.fill_between(x, low + shift, high + shift, color=colour, alpha=0.12, lw=0)
    ax.plot(x, line, color=colour, lw=1.6, marker="o", markersize=2.5, label=label)


# Signatures whose monostatic value is set by the target, not by a numerical cross-pol floor;
# only their change from the monostatic value is meaningful
MONOSTATIC_DEFINED = {
    "span_db",
    "copol_ratio_db",
    "copol_phase_difference_deg",
    ("pauli_fractions_db", 0),
    "entropy",
    "mean_alpha_deg",
}


def figure(pooled, centre, transform, title, path, panels=PANELS):
    columns = 4
    rows = -(-len(panels) // columns)
    fig, axes = plt.subplots(
        rows, columns, figsize=(22, 4.2 * rows), constrained_layout=True, squeeze=False
    )
    for ax, (panel, unit, key, circular, domain) in zip(axes.ravel(), panels):
        for colour, (stem, (features, angles)) in zip(COLOURS, pooled.items()):
            values = transform(values_of(features, key), circular)
            draw(ax, angles, values, centre, colour, label_of(stem), circular)
        style_axis(ax, panel, unit)
        ax.set_xlabel("bistatic angle (deg)", fontsize=8.5, color=MUTED)
        if circular:
            ax.set_ylim(-180, 180)
            ax.set_yticks(np.arange(-180, 181, 90))
        elif domain and transform is identity:
            ax.set_ylim(*domain)
    for ax in axes.ravel()[len(panels) :]:
        ax.set_visible(False)
    axes.ravel()[0].legend(fontsize=7.5, frameon=False)
    fig.suptitle(title, fontsize=10, color=INK, x=0.01, ha="left")
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"Wrote {path}")


def identity(values, circular):
    return values


def change_from_monostatic(values, circular):
    """Signature minus its value at the smallest bistatic angle (the monostatic receiver)."""
    change = values - values[:1]
    return wrap_deg(change) if circular else change


def write_csv(pooled, centre, path):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["design", "bistatic_angle_deg", "panel", "f0", "band_min", "band_max"])
        for stem, (features, angles) in pooled.items():
            for panel, _, key, circular, _ in PANELS:
                values = values_of(features, key)
                if values.ndim == 1:
                    rows = zip(angles, values, values, values)
                else:
                    low, high = band_limits(values, centre, circular)
                    rows = zip(angles, values[:, centre], low, high)
                for angle, f0, low_value, high_value in rows:
                    writer.writerow([label_of(stem), angle, panel, f0, low_value, high_value])
    print(f"Wrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--aim", type=int, default=1, choices=(0, 1))
    parser.add_argument("--out", default=str(SANDBOX_DIR / "summary"))
    parser.add_argument("--designs", nargs="*", default=["MultiRx"])
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    stems = sorted(
        p.stem
        for p in EXPORTS.glob("*.npz")
        if f"rxAim{args.aim}" in p.stem and any(p.stem.startswith(d) for d in args.designs)
    )
    pooled = load(stems)
    some = np.load(EXPORTS / f"{stems[0]}.npz")["frequencies_hz"]
    centre = len(some) // 2
    band = (
        f"line {some[centre] / 1e9:.2f} GHz, band min–max over "
        f"{some[0] / 1e9:.0f}–{some[-1] / 1e9:.0f} GHz"
    )
    aim = "aimed at the target" if args.aim else "on boresight"
    figure(
        pooled,
        centre,
        identity,
        f"Signatures against the bistatic angle, receivers {aim} ({band})",
        out / f"summary_signatures_aim{args.aim}.png",
    )
    figure(
        pooled,
        centre,
        change_from_monostatic,
        f"Change from the monostatic value, receivers {aim} ({band}; phases wrapped; only "
        "signatures whose monostatic value is not a numerical cross-pol floor)",
        out / f"summary_change_aim{args.aim}.png",
        [p for p in PANELS if p[2] in MONOSTATIC_DEFINED],
    )
    write_csv(pooled, centre, out / f"summary_aim{args.aim}.csv")

    pairs = [
        (s, s.replace("rxAim1", "rxAim0"))
        for s in stems
        if args.aim == 1 and (EXPORTS / f"{s.replace('rxAim1', 'rxAim0')}.npz").exists()
    ]
    if pairs:
        difference = {}
        for aimed, boresight in pairs:
            f1, a1 = load([aimed])[aimed]
            f0, a0 = load([boresight])[boresight]
            if not np.allclose(a1, a0):
                continue
            diff = {}
            for _, _, key, circular, _ in PANELS:
                name = key if isinstance(key, str) else key[0]
                if name in diff:
                    continue
                d = f1[name] - f0[name]
                diff[name] = wrap_deg(d) if "phase" in name else d
            difference[aimed] = (diff, a1)
        figure(
            difference,
            centre,
            identity,
            f"rxAim 1 minus rxAim 0: receivers aimed at the target against on boresight ({band})",
            out / "summary_aim_difference.png",
        )


if __name__ == "__main__":
    main()
