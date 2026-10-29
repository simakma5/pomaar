#!/usr/bin/env python3
"""
Architecture A (single-pol, power-only) on the existing exports (host side, NumPy only): the
bistatic leakage into the power-only Pauli powers of single-polarized elements in four linear
states (H, V, +-45 deg), against the coherent Pauli powers of architecture B (dual-pol:
co-located dual-polarized pairs):

    python src/pomaar/simulator/investigation/single_pol_power_only.py [--out <folder>]

With ideal parametric beams a +-45 deg element at a site is (H +- V) / sqrt(2) with the same
beam, so its voltages follow exactly from the H/V matrix S of the same Tx-Rx geometry:
V_qp = p_q^T S p_p (`diagonal_basis_matrix`). Every geometry is first calibrated in power on the
sphere of the same scene at the same geometry (the scalar amplitude calibration each single-pol
pair-array needs: it removes the beam gains and the range spread between receivers and leaves the
target). Everything is compared as Pauli fractions of the span, the features to classify with, so
that the aperture lobe's common drop of all powers, which every architecture sees, does not count
as a polarimetric change: f = |x|^2 / span (x = a, b, c; span = sum of the four H/V powers =
|a|^2 + |b|^2 + |c|^2 + |d|^2) from the complex S, and f^ = the power-only estimate of |x|^2
(`power_only_pauli_powers`, from one cross-pol channel of each basis, the proposal's formulas,
or both) over
the span measured from the four H/V powers (as `pauli_power_fractions` does for B; a sum of
powers, which cannot cancel, unlike the sum of the estimates). With all powers from one pair and
both cross-pol channels the two spans are equal. Differences of fractions are in dB
(-30 dB = 0.001).

Three-receiver scenes (`MultiRx*_rxAim1`: one Tx, receivers at 0, beta and 2 beta), pooled over
the pair bistatic angle like `summarize.py`, per Pauli component:
* the coherent fraction f(beta) of the pair (what B measures there);
* |f(beta) - f(0)|: the target's own change from monostatic, which every architecture sees;
* |f^(beta) - f(beta)| with all eight powers from the same pair, from one cross-pol channel of
  each basis and from both: what the power-only formulas add (both cross-pol channels:
  |d|^2 / span, the antisymmetric power).

All of A's powers come from the same Tx-Rx pair as B's: this is a bistatic-angle analysis and
assumes no layout. Neither the split-pair scenes (a specific interleaved layout, `pairs.py`) nor
any assignment of A's powers to different pairs is used.

Outputs (in <results>/single_pol_power_only/): `three_receivers_<design>.png` per target design
(a row per roll, a column per quantity, a colour per Pauli component), `three_receivers.csv`
(centre frequency and band limits) and `at_angles.csv` (a few angles).
Figures are target-agnostic: no expected values or thresholds; values below -120 dB are drawn
at -120 dB.
"""

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.simulator.array_templates import polarimetric_quads  # noqa: E402
from pomaar.simulator.bistaticity_analysis import (  # noqa: E402
    INK,
    MUTED,
    band_limits,
    receiver_bistatic_angles,
    style_axis,
)
from pomaar.simulator.investigation.project import RESULTS_DIR  # noqa: E402
from pomaar.simulator.investigation.summarize import (  # noqa: E402
    COLOURS,
    EXPORTS,
    draw,
)
from pomaar.simulator.polarimetry_processor import (  # noqa: E402
    bistatic_path,
    channel_ports,
    db10,
    diagonal_basis_matrix,
    pauli_vector,
    power_only_pauli_powers,
    to_scattering_matrices,
)

COMPONENTS = ("a", "b", "c")
COMPONENT_NAMES = {
    "a": r"$|a|^2$ (odd bounce)",
    "b": r"$|b|^2$ (even bounce)",
    "c": r"$|c|^2$ (cross-pol)",
}
FLOOR_DB = -120.0
REPORT_ANGLES_DEG = (1.0, 2.0, 4.0)
# The targets, (design, roll in deg or None). The dihedral at roll 67.5 and 90 deg gives the same
# fractions as at 22.5 and 0 deg (H and V swapped), with less leakage
TARGETS = (
    ("Sphere", None),
    ("Plate", None),
    ("Dihedral", 0.0),
    ("Dihedral", 22.5),
    ("Dihedral", 45.0),
    ("Trihedral", 0.0),
    ("Trihedral", 45.0),
    ("Trihedral", 90.0),
)
# --- loading and calibration --------------------------------------------------------------------


def receivers(path):
    """
    S[rx, tx] of every receiver of an export, sorted by bistatic angle: a list of
    (matrices (n_var, n_freq, 2, 2), bistatic angles (n_var,), two-way paths (n_var,)), plus
    the frequencies.
    """
    data = np.load(path)
    expressions = data["expressions"]
    metadata = json.loads(str(data["metadata"]))
    quads = polarimetric_quads([port for pair in channel_ports(expressions) for port in pair])
    out = [
        (
            to_scattering_matrices(data["s"], expressions, quad),
            receiver_bistatic_angles(metadata, ports),
            port_path(metadata, ports["RxH"], ports["TxH"]),
        )
        for quad, ports in quads.items()
    ]
    out.sort(key=lambda r: r[1].mean())
    return out, data["frequencies_hz"]


def port_path(metadata, rx_port, tx_port):
    """Two-way path (m) Tx -> target -> Rx per variation, from the export's scene geometry."""
    positions = metadata["positions_m"]
    return bistatic_path(
        np.asarray(positions[tx_port]), np.asarray(positions[rx_port]), metadata["target_m"]
    )


def sphere_gain(sphere_matrices, monostatic):
    """Co-pol power of the sphere at each geometry relative to its monostatic one (per freq)."""

    def copol(s):
        return (np.abs(s[..., 0, 0]) ** 2 + np.abs(s[..., 1, 1]) ** 2) / 2

    return copol(sphere_matrices) / copol(monostatic)


def calibrated(matrices, gain):
    return matrices / np.sqrt(gain)[..., None, None]


def coherent(s):
    """|a|^2, |b|^2, |c|^2 and |d|^2 of the Pauli vector, (..., 4)."""
    return np.abs(pauli_vector(s)) ** 2


def coherent_fractions(s):
    """Pauli fractions f (..., 3): |a|^2, |b|^2, |c|^2 over the span (|d|^2 is the rest)."""
    powers = coherent(s)
    return powers[..., :3] / powers.sum(axis=-1, keepdims=True)


def estimated_fractions(s, cross_pol):
    """
    Power-only fractions f^ (..., 3): the estimated |a|^2, |b|^2, |c|^2 over the span measured
    from the H/V powers, which (unlike the sum of the estimates) is a sum of powers and cannot
    cancel.
    """
    hv = np.abs(s) ** 2
    estimate = power_only_pauli_powers(hv, np.abs(diagonal_basis_matrix(s)) ** 2, cross_pol)
    return estimate / hv.sum(axis=(-2, -1))[..., None]


def to_db(fraction):
    """A fraction or a difference of fractions in dB, floored at FLOOR_DB."""
    return np.maximum(db10(fraction), FLOOR_DB)


# --- three-receiver scenes ----------------------------------------------------------------------


def three_receivers(stem, sphere):
    """Products of one MultiRx export against the bistatic angle; see the module docstring."""
    rx, frequencies = receivers(EXPORTS / f"{stem}.npz")
    if len(rx) != 3:
        raise ValueError(f"{stem}: expected three receivers, found {len(rx)}")
    mono_sphere = sphere[0][0][0]
    s = [calibrated(r[0], sphere_gain(sr[0], mono_sphere)) for r, sr in zip(rx, sphere)]
    angles = [r[1] for r in rx]
    if not all(np.allclose(r[1], sr[1], atol=1e-6) for r, sr in zip(rx, sphere)):
        raise ValueError(f"{stem}: the sphere's receivers sit at other bistatic angles")
    reference = coherent_fractions(s[0][0])

    # single pairs, pooled over the angle (the monostatic receiver once)
    pooled = np.concatenate([s[0][:1], s[1], s[2]])
    beta = np.concatenate([angles[0][:1], angles[1], angles[2]])
    order = np.argsort(beta, kind="stable")
    keep = [i for n, i in enumerate(order) if n == 0 or beta[i] - beta[order[n - 1]] > 1e-6]
    pooled, beta = pooled[keep], beta[keep]
    f = coherent_fractions(pooled)
    out = {
        "beta": beta,
        "coherent": to_db(f),
        "change": to_db(np.abs(f - reference)),
        "leak_one_cross": to_db(np.abs(estimated_fractions(pooled, "one") - f)),
        "leak_both_cross": to_db(np.abs(estimated_fractions(pooled, "both") - f)),
        "frequencies": frequencies,
    }
    return out


# --- figures and tables -------------------------------------------------------------------------


# Columns (key, title): f is the coherent Pauli fraction (B), f^ the power-only estimate (A)
THREE_RECEIVER_COLUMNS = (
    ("coherent", "B (coherent): fraction $f$ of the pair"),
    ("change", "Target's own change: $|f(\\beta) - f(0)|$\n(seen by every architecture)"),
    ("leak_one_cross", "A, powers from the same pair, one cross-pol channel per basis\n"
                       "(the proposal's formulas): $|\\hat f - f|$"),
    ("leak_both_cross", "A, powers from the same pair, both cross-pol channels per basis:\n"
                        "$|\\hat f - f|$"),
)  # fmt: skip
# The key to the symbols, printed above both figures
DEFINITIONS = (
    r"$a = (S_\mathrm{HH} + S_\mathrm{VV})/\sqrt{2}$ (odd bounce), "
    r"$b = (S_\mathrm{HH} - S_\mathrm{VV})/\sqrt{2}$ (even bounce), "
    r"$c = (S_\mathrm{HV} + S_\mathrm{VH})/\sqrt{2}$ (cross-pol), "
    r"$d = (S_\mathrm{HV} - S_\mathrm{VH})/\sqrt{2}$ (antisymmetric, 0 monostatically): Pauli "
    r"components of a pair's scattering matrix $S$",
    r"$f = |x|^2 / \mathrm{span}$ for $x = a, b, c$, with span $= |S_\mathrm{HH}|^2 + "
    r"|S_\mathrm{HV}|^2 + |S_\mathrm{VH}|^2 + |S_\mathrm{VV}|^2 = |a|^2 + |b|^2 + |c|^2 + |d|^2$: "
    r"Pauli fraction from the complex $S$ of one Tx–Rx pair, as architecture B (dual-pol, "
    r"coherent) measures it",
    r"$\hat f$: the same fraction as architecture A (single-pol, power-only) estimates it from "
    r"eight powers: $|x|^2$ from $|S_\mathrm{HV}|^2, |S_\mathrm{VH}|^2$ and the $\pm 45^\circ$ "
    r"powers $|V_{++}|^2, |V_{--}|^2, |V_{-+}|^2, |V_{+-}|^2$, over the span from the four H/V "
    r"powers",
    r"    both cross-pol channels per basis: $S_\mathrm{HV}$ and $S_\mathrm{VH}$, $V_{-+}$ and "
    r"$V_{+-}$; one cross-pol channel per basis: the proposal's formulas, $|V_{-+}|^2$ and "
    r"$|S_\mathrm{HV}|^2$ only",
    r"dB $= 10 \log_{10}$ of a fraction or of a difference of fractions: 0 dB = 1 (all of the "
    r"span, or an error as large as all of it), $-10$ dB = 0.1, $-30$ dB = 0.001",
    r"$\beta$: bistatic angle Tx–target–Rx of a pair. Data: ideal identical beams with co-located "
    r"H and V, so the $\pm 45^\circ$ voltages $V_{qp} = \mathbf{p}_q^\mathsf{T} S\, \mathbf{p}_p$ "
    r"follow exactly from $S$; every pair is calibrated on the sphere at its own geometry",
)


def band_text(frequencies, centre):
    return (
        f"Lines {frequencies[centre] / 1e9:.2f} GHz, bands min–max over "
        f"{frequencies[0] / 1e9:.0f}–{frequencies[-1] / 1e9:.0f} GHz, values below "
        f"{FLOOR_DB:.0f} dB drawn at {FLOOR_DB:.0f} dB"
    )


def header(subfigure, title, sketch, frequencies):
    """Scene sketch (left) and the key to the symbols (right) in the figure's top band."""
    sketch_ax, text_ax = subfigure.subplots(1, 2, width_ratios=(1, 3.4))
    sketch(sketch_ax)
    handles = [
        matplotlib.lines.Line2D([], [], color=colour, lw=1.6, marker="o", markersize=2.5,
                                label=COMPONENT_NAMES[component])
        for colour, component in zip(COLOURS, COMPONENTS)
    ]  # fmt: skip
    text_ax.legend(handles=handles, loc="lower left", ncols=3, fontsize=9, frameon=False,
                   title="Colours (every panel):", title_fontsize=9, alignment="left")  # fmt: skip
    text_ax.axis("off")
    lines = DEFINITIONS + (band_text(frequencies, len(frequencies) // 2),)
    text_ax.text(0, 1, "\n".join(lines), va="top", ha="left", fontsize=9, color=INK,
                 linespacing=1.7, transform=text_ax.transAxes)  # fmt: skip
    subfigure.suptitle(title, fontsize=11, color=INK, x=0.005, ha="left")


def antenna(ax, x, label, colour=INK):
    ax.plot(x, 0, marker="^", markersize=8, color=colour, linestyle="none")
    ax.annotate(label, (x, 0), xytext=(0, -11), textcoords="offset points", ha="center",
                va="top", fontsize=8, color=colour)  # fmt: skip


def scene_axes(ax, title):
    ax.plot(0, 1, marker="s", markersize=9, color=INK)
    ax.annotate("target, 4 m on boresight", (0, 1), xytext=(8, 0), textcoords="offset points",
                va="center", fontsize=8, color=INK)  # fmt: skip
    ax.set_title(title, loc="left", fontsize=9, color=INK)
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-0.45, 1.12)
    ax.axis("off")


def sketch_three_receivers(ax):
    xs = {"Tx, Rx$_0$\n(H, V)": 0.0, "Rx$_1$ (H, V)\nangle $\\beta$": 0.4,
          "Rx$_2$ (H, V)\nangle $2\\beta$": 0.85}  # fmt: skip
    for x in xs.values():
        ax.plot([x, 0], [0, 1], color=MUTED, lw=1)
    for label, x in xs.items():
        antenna(ax, x, label)
    scene_axes(ax, "Scene (not to scale): one Tx, three receivers\naimed at the target; "
                   "each pair is analysed on its own")  # fmt: skip


def pauli_figure(results, columns, title, sketch, path):
    """A row per target, a column per quantity, a colour per Pauli component."""
    frequencies = next(iter(results.values()))["frequencies"]
    centre = len(frequencies) // 2
    fig = plt.figure(figsize=(5 * len(columns) + 2, 3 + 3.4 * len(results)), layout="constrained")
    top, bottom = fig.subfigures(2, 1, height_ratios=(3, 3.4 * len(results)))
    header(top, title, sketch, frequencies)
    axes = bottom.subplots(len(results), len(columns), squeeze=False, sharey="row")
    for i, (row, (name, r)) in enumerate(zip(axes, results.items())):
        for j, (ax, (key, column_title)) in enumerate(zip(row, columns)):
            for colour, component in zip(COLOURS, range(len(COMPONENTS))):
                label = COMPONENT_NAMES[COMPONENTS[component]]
                draw(ax, r["beta"], r[key][..., component], centre, colour, label, False)
            style_axis(ax, column_title if i == 0 else "", f"{name}\ndB" if j == 0 else "")
            ax.set_xlabel("pair angle $\\beta$ (deg)", fontsize=8.5, color=MUTED)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    print(f"Wrote {path}")


def stem_of(prefix, scene, target):
    """('MultiRx', 'rxAim1', ('Dihedral', 22.5)) -> 'MultiRxDihedral_rxAim1_targetRoll22.5'."""
    design, roll = target
    return f"{prefix}{design}_{scene}" + ("" if roll is None else f"_targetRoll{roll:g}")


def label_of(target):
    design, roll = target
    return design if roll is None else f"{design}, roll {roll:g}°"


def three_receiver_figure(results, path):
    named = {label_of(target): r for target, r in results.items()}
    title = (
        "Architecture A (single-pol, power-only) against B (dual-pol, coherent): three-receiver "
        "scenes"
    )
    pauli_figure(named, THREE_RECEIVER_COLUMNS, title, sketch_three_receivers, path)


def write_long_csv(results, keys, path):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target", "roll_deg", "quantity", "component", "beta_deg", "f0_db",
                         "band_min_db", "band_max_db"])  # fmt: skip
        for target, r in results.items():
            centre = len(r["frequencies"]) // 2
            for key in keys:
                for component, name in enumerate(COMPONENTS):
                    values = r[key][..., component]
                    low, high = band_limits(values, centre)
                    for x, f0, lo, hi in zip(r["beta"], values[:, centre], low, high):
                        writer.writerow([*target, key, name, x, f0, lo, hi])
    print(f"Wrote {path}")


def at_angles(rows, results, keys):
    for target, r in results.items():
        centre = len(r["frequencies"]) // 2
        for angle in REPORT_ANGLES_DEG:
            for key in keys:
                x = r["beta"]
                i = int(np.argmin(np.abs(x - angle)))
                if abs(x[i] - angle) > 1e-3:
                    continue
                for component, name in enumerate(COMPONENTS):
                    values = r[key][i, :, component]
                    rows.append([*target, angle, key, name, values[centre],
                                 values.max()])  # fmt: skip


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default=str(RESULTS_DIR / "single_pol_power_only"))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    sphere, _ = receivers(EXPORTS / "MultiRxSphere_rxAim1.npz")
    ladder = {target: three_receivers(stem_of("MultiRx", "rxAim1", target), sphere)
              for target in TARGETS}  # fmt: skip
    for design in dict.fromkeys(design for design, _ in TARGETS):
        family = design.lower()
        three_receiver_figure({t: r for t, r in ladder.items() if t[0] == design},
                              out / f"three_receivers_{family}.png")  # fmt: skip
    keys = [key for key, _ in THREE_RECEIVER_COLUMNS]
    write_long_csv(ladder, keys, out / "three_receivers.csv")
    rows = []
    at_angles(rows, ladder, keys)
    with open(out / "at_angles.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target", "roll_deg", "beta_deg", "quantity", "component",
                         "f0_db", "band_max_db"])  # fmt: skip
        writer.writerows(rows)
    print(f"Wrote {out / 'at_angles.csv'}")


if __name__ == "__main__":
    main()
