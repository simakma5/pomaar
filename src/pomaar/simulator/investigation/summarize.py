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
* `summary_aim<a>.csv`: the centre-frequency values and band limits behind the figures;
* `summary_departures_aim<a>.csv`: per design, the first bistatic angle at which each
  descriptor departs from its monostatic value by the levels of `DEPARTURES`, at the centre
  frequency and anywhere in the band; cross-pol descriptors count only where both cross-pol
  channels stand clear of the numerical floor (`CROSS_POL_DEFINED_DB`). The levels are an
  interpretation (how much change counts) and appear only in this table, never in the figures.
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


# (descriptor, feature key, kind, levels): 'change' is |value - monostatic value| (wrapped for
# phases), 'below' a value that falls under the level, 'above' one that rises over it
DEPARTURES = [
    ("span (lobe)", "span_db", "drop", (3.0, 10.0)),
    ("co-pol ratio |VV/HH| (dB)", "copol_ratio_db", "change", (0.5, 1.0)),
    ("co-pol phase difference (deg)", "copol_phase_difference_deg", "change", (5.0, 10.0)),
    ("co/cross Tx H |HH/VH| (dB)", "co_minus_cross_txh_db", "below", (30.0, 20.0)),
    ("co/cross Tx V |VV/HV| (dB)", "co_minus_cross_txv_db", "below", (30.0, 20.0)),
    ("cross-pol ratio |HV/VH| (dB)", "reciprocity_db", "absolute", (1.0, 3.0)),
    ("entropy H", "entropy", "above", (0.1, 0.3)),
    ("mean alpha (deg)", "mean_alpha_deg", "above", (10.0, 20.0)),
]


# Cross-pol descriptors are only defined where both cross-pol channels stand clear of the
# numerical cross-pol floor (about 60 dB below co-pol at 4 rays per wavelength, trust check T4)
CROSS_POL_DEFINED_DB = 50.0
CROSS_POL_KEYS = {"reciprocity_db", "reciprocity_phase_deg", "asymmetry_db"}


def cross_pol_defined(features, index=None):
    """Mask of the angles where both co/cross ratios are below CROSS_POL_DEFINED_DB."""
    h, v = features["co_minus_cross_txh_db"], features["co_minus_cross_txv_db"]
    if index is not None:
        h, v = h[:, index], v[:, index]
    return (h < CROSS_POL_DEFINED_DB) & (v < CROSS_POL_DEFINED_DB)


def first_departure(angles, values, kind, level, circular, defined=None):
    """
    First bistatic angle at which `values` (per angle) departs by `level`, among the angles
    where `defined` is true (all by default); None if never.
    """
    if kind == "change":
        deviation = values - values[0]
        exceeded = np.abs(wrap_deg(deviation) if circular else deviation) > level
    elif kind == "drop":
        exceeded = values[0] - values > level
    elif kind == "below":
        exceeded = values < level
    elif kind == "absolute":
        exceeded = np.abs(values) > level
    else:
        exceeded = values > level
    exceeded[0] = False
    if defined is not None:
        exceeded &= defined
    return round(float(angles[np.argmax(exceeded)]), 2) if exceeded.any() else None


def write_departures(pooled, centre, path):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["design", "descriptor", "level", "first_angle_f0_deg", "first_angle_band_deg"]
        )
        for stem, (features, angles) in pooled.items():
            for descriptor, key, kind, levels in DEPARTURES:
                values = features[key]
                circular = "phase" in key
                cross = key in CROSS_POL_KEYS

                def mask(index):
                    return cross_pol_defined(features, index) if cross else None

                for level in levels:
                    if values.ndim == 1:
                        at_f0 = in_band = first_departure(angles, values, kind, level, circular)
                    else:
                        at_f0 = first_departure(
                            angles, values[:, centre], kind, level, circular, mask(centre)
                        )
                        # anywhere in the band: the earliest angle over all frequencies
                        firsts = [
                            first_departure(angles, values[:, i], kind, level, circular, mask(i))
                            for i in range(values.shape[1])
                        ]
                        firsts = [f for f in firsts if f is not None]
                        in_band = min(firsts) if firsts else None
                    writer.writerow([label_of(stem), descriptor, level, at_f0, in_band])
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
    write_departures(pooled, centre, out / f"summary_departures_aim{args.aim}.csv")

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
