#!/usr/bin/env python3
"""
Differential reflectivity from split antenna pairs (host side, NumPy only), from the exports of
`build_scene --layout pairs` (one per target and gapRatio, swept over pairAngle):

    python src/pomaar/simulator/investigation/pairs.py [--designs PairsSphere ...] [--out <folder>]

VV is measured by the pair TxV-RxV and HH by the pair TxH-RxH, both of base b, with their
bisectors (b + g) / R apart; the co-located dual-polarized pair of the same base (bisector on
boresight) gives the reference. Per pair bistatic angle and frequency:

* Z_DR = |S_HH|^2 / |S_VV|^2 (dB) and the co-pol phase difference arg(S_HH S_VV*) of the split
  pairs and of the reference pair;
* their differences, split minus reference: what measuring HH and VV with separate pairs adds.

Outputs (in <sandbox>/pairs/): `pairs_zdr.png` and `pairs_zdr.csv` (centre frequency, band
limits). The figure is target-agnostic: no expected values or thresholds.
"""

import argparse
import csv
import json
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.simulator.array_templates import parse_port  # noqa: E402
from pomaar.simulator.bistaticity_analysis import INK, style_axis, wrap_deg  # noqa: E402
from pomaar.simulator.investigation.sandbox import SANDBOX_DIR  # noqa: E402
from pomaar.simulator.investigation.summarize import COLOURS, EXPORTS, draw  # noqa: E402
from pomaar.simulator.polarimetry_processor import (  # noqa: E402
    bistatic_angle_deg,
    channel_ports,
    db20,
)

# (role, pair key, polarization) of the channels used
SPLIT_VV = (("Rx", "A", "V"), ("Tx", "A", "V"))
SPLIT_HH = (("Rx", "B", "H"), ("Tx", "B", "H"))
REF_VV = (("Rx", "R", "V"), ("Tx", "R", "V"))
REF_HH = (("Rx", "R", "H"), ("Tx", "R", "H"))


def port_id(port):
    role, pol, key = parse_port(port)
    return role, key.strip("_").split("_")[1], pol


def channel(data, ports):
    """S of the channel (Rx, Tx) given by `ports` ((role, key, pol) each): (n_var, n_freq)."""
    pairs = [tuple(port_id(p) for p in pair) for pair in channel_ports(data["expressions"])]
    return data["s"][:, pairs.index(ports)]


def load(path):
    """Per export: pair bistatic angle, bisector separation, Z_DR and phases of split and ref."""
    data = dict(np.load(path))
    metadata = json.loads(str(data["metadata"]))
    positions = {port_id(p): np.asarray(v) for p, v in metadata["positions_m"].items()}
    target = np.asarray(metadata["target_m"])
    tx_v, rx_v = positions[SPLIT_VV[1]], positions[SPLIT_VV[0]]
    tx_h, rx_h = positions[SPLIT_HH[1]], positions[SPLIT_HH[0]]

    def bisector(tx, rx):
        u = (tx - target) / np.linalg.norm(tx - target, axis=-1, keepdims=True)
        v = (rx - target) / np.linalg.norm(rx - target, axis=-1, keepdims=True)
        return (u + v) / np.linalg.norm(u + v, axis=-1, keepdims=True)

    separation = np.degrees(
        np.arccos(np.clip(np.sum(bisector(tx_v, rx_v) * bisector(tx_h, rx_h), axis=-1), -1, 1))
    )
    angle_v, angle_h = (
        bistatic_angle_deg(tx_v, rx_v, target),
        bistatic_angle_deg(tx_h, rx_h, target),
    )
    s = {name: channel(data, ports) for name, ports in (
        ("split_vv", SPLIT_VV), ("split_hh", SPLIT_HH), ("ref_vv", REF_VV), ("ref_hh", REF_HH)
    )}  # fmt: skip
    out = {
        "pair_angle_deg": angle_v,
        "pair_angle_mismatch_deg": np.abs(angle_v - angle_h),
        "bisector_separation_deg": separation,
        "gap_ratio": float(metadata["nominal_variables"]["gapRatio"]),
        "frequencies_hz": data["frequencies_hz"],
    }
    for kind in ("split", "ref"):
        hh, vv = s[f"{kind}_hh"], s[f"{kind}_vv"]
        out[f"zdr_{kind}_db"] = db20(hh) - db20(vv)
        out[f"hh_vv_phase_{kind}_deg"] = np.degrees(np.angle(hh * np.conj(vv)))
    out["zdr_split_minus_ref_db"] = out["zdr_split_db"] - out["zdr_ref_db"]
    out["hh_vv_phase_split_minus_ref_deg"] = wrap_deg(
        out["hh_vv_phase_split_deg"] - out["hh_vv_phase_ref_deg"]
    )
    order = np.argsort(angle_v)
    return {k: v[order] if np.ndim(v) and len(v) == len(order) else v for k, v in out.items()}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--designs", nargs="*", default=None, help="Design stems (default: all)")
    parser.add_argument("--out", default=str(SANDBOX_DIR / "pairs"))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    paths = sorted(EXPORTS.glob("Pairs*.npz"))
    # one series per gapRatio: the design is the export stem without its gapRatio token
    design_of = {p: re.sub(r"_gapRatio[\d.]+", "", p.stem) for p in paths}
    designs = sorted(set(design_of.values()))
    if args.designs:
        designs = [d for d in designs if d.split("_")[0] in args.designs or d in args.designs]
    results = {d: [load(p) for p in paths if design_of[p] == d] for d in designs}
    for runs in results.values():
        runs.sort(key=lambda r: r["gap_ratio"])
    centre = len(next(iter(results.values()))[0]["frequencies_hz"]) // 2
    frequencies = next(iter(results.values()))[0]["frequencies_hz"]

    panels = [
        (
            r"$Z_\mathrm{DR} = |S_\mathrm{HH}|^2 / |S_\mathrm{VV}|^2$, split pairs",
            "dB",
            "zdr_split_db",
        ),
        (r"$Z_\mathrm{DR}$, co-located reference pair", "dB", "zdr_ref_db"),
        (r"$\Delta Z_\mathrm{DR}$: split minus reference", "dB", "zdr_split_minus_ref_db"),
        (
            r"$\Delta\angle(S_\mathrm{HH} S_\mathrm{VV}^*)$: split minus reference",
            "deg",
            "hh_vv_phase_split_minus_ref_deg",
        ),
    ]
    fig, axes = plt.subplots(
        len(designs), len(panels), figsize=(5.2 * len(panels), 3.6 * len(designs)),
        constrained_layout=True, squeeze=False,
    )  # fmt: skip
    for row, design in zip(axes, designs):
        for ax, (title, unit, key) in zip(row, panels):
            for colour, run in zip(COLOURS, results[design]):
                separation = run["bisector_separation_deg"] / np.maximum(
                    run["pair_angle_deg"], 1e-12
                )
                apart = np.median(separation)
                label = rf"$g = {run['gap_ratio']:g}\,b$ (bisectors {apart:.1f} $\beta$ apart)"
                draw(ax, run["pair_angle_deg"], run[key], centre, colour, label, "phase" in key)
            name = design.removeprefix("Pairs").replace("_targetRoll", ", roll ")
            style_axis(ax, f"{name}: {title}", unit)
            ax.set_xlabel(r"pair bistatic angle $\beta$ (deg)")
        row[0].legend(fontsize=7, frameon=False)
    fig.suptitle(
        "Split antenna pairs (TxV–RxV, gap, TxH–RxH) against a co-located dual-polarized pair of "
        f"the same base. Line {frequencies[centre] / 1e9:.2f} GHz, band min–max over "
        f"{frequencies[0] / 1e9:.0f}–{frequencies[-1] / 1e9:.0f} GHz",
        fontsize=10,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.savefig(out / "pairs_zdr.png", dpi=150)
    plt.close(fig)

    keys = [p[2] for p in panels] + ["hh_vv_phase_split_deg", "hh_vv_phase_ref_deg"]
    with open(out / "pairs_zdr.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["design", "gap_ratio", "pair_angle_deg", "bisector_separation_deg"]
            + [f"{k}@f0" for k in keys]
            + [f"{k}_band_{e}" for k in keys for e in ("min", "max")]
        )
        for design, runs in results.items():
            for run in runs:
                for i, angle in enumerate(run["pair_angle_deg"]):
                    deviation = {
                        k: wrap_deg(run[k][i] - run[k][i, centre]) if "phase" in k
                        else run[k][i] - run[k][i, centre]
                        for k in keys
                    }  # fmt: skip
                    writer.writerow(
                        [design, run["gap_ratio"], angle, run["bisector_separation_deg"][i]]
                        + [run[k][i, centre] for k in keys]
                        + [
                            run[k][i, centre] + f(deviation[k])
                            for k in keys
                            for f in (np.min, np.max)
                        ]
                    )
    mismatch = max(r["pair_angle_mismatch_deg"].max() for runs in results.values() for r in runs)
    print(
        f"Wrote {out / 'pairs_zdr.png'} and the CSV; largest HH/VV pair-angle mismatch "
        f"{mismatch:.2g} deg"
    )


if __name__ == "__main__":
    main()
