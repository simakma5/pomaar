#!/usr/bin/env python3
"""
Reference check of receiver aiming (host side, NumPy only): for every design solved with both
rxAim 0 (receivers on boresight) and rxAim 1 (receivers turned towards the target), the target
matrix of each receiver and bistatic angle is fitted as one complex scalar times its on-boresight
counterpart, S_aim = g S_boresight (least squares over the four channels, per frequency):

    python src/pomaar/simulator/investigation/aim_check.py [--out <folder>]

Every descriptor of `CanonicalTargetProcessor.signatures` except the span is a ratio, a phase
difference or a span-normalized quantity, i.e. unchanged by S -> g S. A small fit residual
||S_aim - g S_boresight|| / ||S_aim|| therefore shows that aiming changes the span by |g|^2 and
nothing else. Outputs (in <results>/rx_aim/):

* `rx_aim_scalar_fit.csv`: per design, receiver and bistatic angle, the residual, |g|^2 and
  arg g at the centre frequency with their band limits, and the span relative to the monostatic
  span (to tell nulls apart);
* `rx_aim_scalar_fit.png`: the same against the bistatic angle, one series per design.
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.simulator.array_templates import polarimetric_quads  # noqa: E402
from pomaar.simulator.bistaticity_analysis import (  # noqa: E402
    INK,
    MUTED,
    receiver_bistatic_angles,
    style_axis,
)
from pomaar.simulator.investigation.project import RESULTS_DIR  # noqa: E402
from pomaar.simulator.investigation.summarize import (  # noqa: E402
    COLOURS,
    EXPORTS,
    draw,
    label_of,
)
from pomaar.simulator.polarimetry_processor import (  # noqa: E402
    CanonicalTargetProcessor,
    channel_ports,
    db10,
    db20,
)


def scalar_fit(boresight, aimed):
    """(g, residual) of aimed = g boresight per matrix (..., 2, 2): least squares over S."""
    g = np.sum(np.conj(boresight) * aimed, axis=(-2, -1)) / np.sum(
        np.abs(boresight) ** 2, axis=(-2, -1)
    )
    residual = np.linalg.norm(aimed - g[..., None, None] * boresight, axis=(-2, -1))
    return g, residual / np.linalg.norm(aimed, axis=(-2, -1))


def fit_design(boresight_path, aimed_path):
    """Per receiver and variation, sorted by bistatic angle: {quantity: array}."""
    expressions = np.load(boresight_path)["expressions"]
    quads = polarimetric_quads([p for pair in channel_ports(expressions) for p in pair])
    parts = {k: [] for k in ("angle", "receiver", "residual", "g", "span")}
    for quad, ports in quads.items():
        boresight = CanonicalTargetProcessor.from_npz(boresight_path, quad=quad)
        aimed = CanonicalTargetProcessor.from_npz(aimed_path, quad=quad)
        angles = receiver_bistatic_angles(aimed.metadata, ports)
        if not np.allclose(angles, receiver_bistatic_angles(boresight.metadata, ports)):
            raise ValueError(f"{aimed_path.name} and {boresight_path.name} differ in geometry.")
        s_boresight, s_aimed = boresight.target_matrices(), aimed.target_matrices()
        g, residual = scalar_fit(s_boresight, s_aimed)
        parts["angle"].append(angles)
        parts["receiver"] += [quad.split("/")[0]] * len(angles)
        parts["residual"].append(residual)
        parts["g"].append(g)
        parts["span"].append(np.sum(np.abs(s_aimed) ** 2, axis=(-2, -1)))
    result = {k: np.concatenate(v) if k != "receiver" else np.array(v) for k, v in parts.items()}
    order = np.argsort(result["angle"], kind="stable")
    result = {k: v[order] for k, v in result.items()}
    monostatic = result["span"][np.argmin(result["angle"])]
    result["span_re_monostatic_db"] = db10(result["span"] / monostatic)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default=str(RESULTS_DIR / "rx_aim"))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    pairs = [
        (p, EXPORTS / p.name.replace("rxAim0", "rxAim1"))
        for p in sorted(EXPORTS.glob("MultiRx*_rxAim0*.npz"))
        if (EXPORTS / p.name.replace("rxAim0", "rxAim1")).exists()
    ]
    fits = {label_of(aimed.stem): fit_design(boresight, aimed) for boresight, aimed in pairs}
    centre = next(iter(fits.values()))["residual"].shape[1] // 2

    with open(out / "rx_aim_scalar_fit.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "design",
                "receiver",
                "bistatic_angle_deg",
                "residual_db@f0",
                "residual_db_band_max",
                "gain_db@f0",
                "gain_db_band_min",
                "gain_db_band_max",
                "gain_phase_deg@f0",
                "gain_phase_deg_band_max_abs",
                "span_re_monostatic_db@f0",
            ]
        )
        for design, fit in fits.items():
            residual, gain = db20(fit["residual"]), db20(np.abs(fit["g"]))
            phase = np.degrees(np.angle(fit["g"]))
            for i, angle in enumerate(fit["angle"]):
                writer.writerow(
                    [
                        design,
                        fit["receiver"][i],
                        angle,
                        residual[i, centre],
                        residual[i].max(),
                        gain[i, centre],
                        gain[i].min(),
                        gain[i].max(),
                        phase[i, centre],
                        np.abs(phase[i]).max(),
                        fit["span_re_monostatic_db"][i, centre],
                    ]
                )

    fig, axes = plt.subplots(1, 4, figsize=(22, 4.6), constrained_layout=True)
    aimed, boresight = r"S_\mathrm{aim}", r"S_\mathrm{boresight}"
    panels = [
        (
            rf"Fit residual $\|{aimed} - g\,{boresight}\|^2 / \|{aimed}\|^2$",
            "dB",
            lambda f: db20(f["residual"]),
        ),
        ("Gain $|g|^2$ of aiming", "dB", lambda f: db20(np.abs(f["g"]))),
        (r"Phase $\arg g$ of aiming", "deg", lambda f: np.degrees(np.angle(f["g"]))),
        (
            r"Span $\|S\|^2$ of the aimed receivers",
            "dB re monostatic",
            lambda f: f["span_re_monostatic_db"],
        ),
    ]
    for ax, (title, unit, quantity) in zip(axes, panels):
        for colour, (design, fit) in zip(COLOURS, fits.items()):
            off = fit["angle"] > 0  # the monostatic receiver is not turned
            values = np.where(np.isfinite(quantity(fit)), quantity(fit), np.nan)
            draw(ax, fit["angle"][off], values[off], centre, colour, design, False)
        style_axis(ax, title, unit)
        ax.set_xlabel("bistatic angle (deg)", fontsize=8.5, color=MUTED)
    axes[0].legend(fontsize=7, frameon=False)
    frequencies = np.load(pairs[0][0])["frequencies_hz"]
    fig.suptitle(
        "Receivers aimed at the target (rxAim 1) against on boresight (rxAim 0): "
        rf"${aimed} = g\,{boresight}$. Line {frequencies[centre] / 1e9:.2f} GHz, band min–max over "
        f"{frequencies[0] / 1e9:.0f}–{frequencies[-1] / 1e9:.0f} GHz",
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(out / "rx_aim_scalar_fit.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {out / 'rx_aim_scalar_fit.csv'} and {out / 'rx_aim_scalar_fit.png'}")


if __name__ == "__main__":
    main()
