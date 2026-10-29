#!/usr/bin/env python3
"""
Size and range scaling of the multi-receiver results (host side, NumPy only): every given
export is pooled by bistatic angle and compared against the bistatic angle beta and against
beta / beta_0, beta_0 being the first null of its span at the centre frequency:

    python src/pomaar/simulator/investigation/scaling.py MultiRxTrihedral_rxAim1_targetRoll0 \
        MultiRxTrihedralS80_rxAim1 ... [--out <folder>]

Outputs (in <results>/scaling/):
* `scaling.png`: span, Z_DR and co-pol phase difference, and each cross-pol channel relative to the
  monostatic span (where it stands EMERGED_DB above its monostatic value, which is the
  numerical floor for these targets), against beta (log axis for the cross-pol) and against
  beta / beta_0;
* `scaling.csv`: per export, beta_0, the far-field distance 2 D^2 / lambda of the target's
  aperture against its range, the power-law fit of the cross-pol over FIT_RANGE_DEG (slope in
  dB per doubling of beta and level at 1 deg) and the departure angles of
  `summarize.DEPARTURES` in deg and in units of beta_0.
"""

import argparse
import ast
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pomaar.simulator.bistaticity_analysis import (  # noqa: E402
    INK,
    LABELS,
    MUTED,
    pool_receivers,
    style_axis,
)
from pomaar.simulator.investigation.project import EXPORTS_DIR, RESULTS_DIR  # noqa: E402
from pomaar.simulator.investigation.summarize import (  # noqa: E402
    COLOURS,
    DEPARTURES,
    defined_mask,
    first_departure,
    monostatic_defined,
)

EXPORTS = EXPORTS_DIR
# the bistatic angles over which the emerging cross-pol is fitted by a power law
FIT_RANGE_DEG = (0.25, 1.5)
EMERGED_DB = 10.0
C0 = 299_792_458.0


def label_of(stem, metadata):
    nominal = metadata["nominal_variables"]
    nominal = ast.literal_eval(nominal) if isinstance(nominal, str) else nominal
    target = next(t for t in ("Trihedral", "Dihedral", "Sphere") if t in stem)
    return (
        f"{target} {nominal['targetSize']}, {nominal['targetRange']}, roll {nominal['targetRoll']}"
    ), nominal


def length_m(text):
    for unit, scale in (("mm", 1e-3), ("cm", 1e-2), ("meter", 1.0), ("m", 1.0)):
        if text.endswith(unit):
            return float(text.removesuffix(unit)) * scale
    return float(text)


def first_null(angles, span_db):
    """First local minimum of the span (deg), refined by a parabola through its neighbours."""
    rising = np.flatnonzero(np.diff(span_db) > 0)
    if not rising.size or rising[0] == 0:
        return None
    i = rising[0]
    x, y = angles[i - 1 : i + 2], span_db[i - 1 : i + 2]
    a, b, _ = np.polyfit(x, y, 2)
    return float(-b / (2 * a)) if a > 0 and x[0] < -b / (2 * a) < x[-1] else float(angles[i])


def power_law(angles, values, defined):
    """Slope (dB per doubling of beta) and level at 1 deg of `values` over FIT_RANGE_DEG."""
    use = (angles >= FIT_RANGE_DEG[0] - 1e-6) & (angles <= FIT_RANGE_DEG[1] + 1e-6) & defined
    if use.sum() < 3:
        return None, None
    slope, level = np.polyfit(np.log2(angles[use]), values[use], 1)
    return float(slope), float(level)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("stems", nargs="+", help="Export names in <exports>")
    parser.add_argument("--out", default=str(RESULTS_DIR / "scaling"))
    args = parser.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows, series = [], []
    for stem in args.stems:
        path = EXPORTS / f"{stem}.npz"
        features, angles, _, _ = pool_receivers(path)
        data = np.load(path)
        metadata = json.loads(str(data["metadata"]))
        label, nominal = label_of(stem, metadata)
        frequencies = data["frequencies_hz"]
        centre = len(frequencies) // 2
        wavelength = C0 / frequencies[centre]
        size, distance = length_m(nominal["targetSize"]), length_m(nominal["targetRange"])
        # aperture across the target as seen from the antenna (trihedral: the edge times
        # sqrt 2, dihedral: twice the face times cos 45 deg; both a sqrt 2)
        far_field = 2 * (size * np.sqrt(2)) ** 2 / wavelength
        span = features["span_db"][:, centre]
        null = first_null(angles, span)
        reference = features["span_db"][0, centre]
        row = {
            "export": stem,
            "label": label,
            "size_mm": size * 1e3,
            "range_m": distance,
            "far_field_m": round(far_field, 2),
            "first_null_deg": None if null is None else round(null, 3),
            "monostatic_span_db": round(float(reference), 2),
        }
        for channel in ("hv", "vh"):
            values = features[f"power_{channel}_db"][:, centre] - reference
            # the emerging part: 10 dB above the channel's monostatic value (its floor here)
            slope, level = power_law(angles, values, values > values[0] + EMERGED_DB)
            row[f"{channel}_slope_db_per_doubling"] = None if slope is None else round(slope, 1)
            row[f"{channel}_at_1deg_db"] = None if level is None else round(level, 1)
        for descriptor, key, kind, levels in DEPARTURES:
            values = features[key]
            circular = "phase" in key
            level = levels[0]
            if kind == "change" and not monostatic_defined(features, key, centre):
                angle = None
            else:
                line = values if values.ndim == 1 else values[:, centre]
                mask = defined_mask(features, key, centre)
                angle = first_departure(angles, line, kind, level, circular, mask)
            row[f"{descriptor} [{level:g}] deg"] = angle
            row[f"{descriptor} [{level:g}] / beta_0"] = (
                None if angle is None or null is None else round(angle / null, 3)
            )
        rows.append(row)
        series.append((label, features, angles, centre, null, reference))
        print(
            f"{label}: beta_0 {row['first_null_deg']} deg, far field {far_field:.1f} m; "
            f"HV {row['hv_slope_db_per_doubling']} dB per doubling, "
            f"{row['hv_at_1deg_db']} dB at 1 deg"
        )

    with open(out / "scaling.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {out / 'scaling.csv'}")

    panels = [
        (f"Span {LABELS['span_db']} minus monostatic", "dB", "span"),
        (LABELS["zdr_db"], "dB", "zdr_db"),
        (LABELS["hh_vv_phase_deg"], "deg", "hh_vv_phase_deg"),
        (f"{LABELS['power_hv_db']} minus monostatic span", "dB", "power_hv_db"),
        (f"{LABELS['power_vh_db']} minus monostatic span", "dB", "power_vh_db"),
    ]
    fig, axes = plt.subplots(2, len(panels), figsize=(5 * len(panels), 8), constrained_layout=True)
    for column, (title, unit, key) in enumerate(panels):
        for colour, (label, features, angles, centre, null, reference) in zip(COLOURS, series):
            if key == "span":
                values = features["span_db"][:, centre] - reference
            elif key.startswith("power_"):
                values = features[key][:, centre] - reference
                values = np.where(values > values[0] + EMERGED_DB, values, np.nan)
            else:
                values = features[key][:, centre]
            for row, x in enumerate((angles, angles / null if null else None)):
                if x is None:
                    continue
                ax = axes[row, column]
                keep = x > 0 if key.startswith("power_") else np.ones_like(x, bool)
                ax.plot(x[keep], values[keep], color=colour, lw=1.4, marker="o", markersize=2.5,
                        label=label)  # fmt: skip
        for row, xlabel in enumerate(
            (r"bistatic angle $\beta$ (deg)", r"$\beta / \beta_0$ (first null)")
        ):
            ax = axes[row, column]
            style_axis(ax, title, unit)
            ax.set_xlabel(xlabel, fontsize=8.5, color=MUTED)
            if key.startswith("power_"):
                ax.set_xscale("log")
            if "phase" in key:
                ax.set_ylim(-180, 180)
                ax.set_yticks(np.arange(-180, 181, 90))
    axes[0, 0].legend(fontsize=7.5, frameon=False)
    fig.suptitle(
        "Size and range scaling at the centre frequency (cross-pol only where 10 dB above its "
        "monostatic value)",
        fontsize=10, color=INK, x=0.01, ha="left",
    )  # fmt: skip
    fig.savefig(out / "scaling.png", dpi=140)
    plt.close(fig)
    print(f"Wrote {out / 'scaling.png'}")


if __name__ == "__main__":
    main()
