#!/usr/bin/env python3
"""
Trust checks of the multi-receiver scene (host side, NumPy only), from the sandbox exports:

    python src/pomaar/simulator/investigation/trust_checks.py level <export> [--range 4]
    python src/pomaar/simulator/investigation/trust_checks.py paths <export> --bistatic-angle 1
    python src/pomaar/simulator/investigation/trust_checks.py crosspol <export> [<export> ...]
    python src/pomaar/simulator/investigation/trust_checks.py strip <export>
    python src/pomaar/simulator/investigation/trust_checks.py labels <dihedral> <trihedral>

`level` compares the monostatic trihedral return with the radar equation, `paths` the two-way
path from the phase slope with the trihedral's geometric-optics path, `crosspol` the cross-pol
floor of every receiver, `strip` the co-pol levels of a PEC strip along H (uninformative:
physical optics ignores the polarization of a flat plate at normal incidence), `labels` the
polarization labels: relative to a trihedral, a dihedral with its seam along H (targetRoll 0)
flips the sign of the component perpendicular to the seam, i.e. of S_VV, and keeps S_HH.
"""

import argparse
from pathlib import Path

import numpy as np

from pomaar.simulator.investigation.sandbox import SANDBOX_DIR
from pomaar.simulator.polarimetry_processor import SPEED_OF_LIGHT, CanonicalTargetProcessor, db20

EXPORTS = SANDBOX_DIR / "exports"


def export_path(name):
    path = Path(name)
    return path if path.suffix == ".npz" else EXPORTS / f"{name}.npz"


def receivers(path):
    """{quad name: matrices (variation x frequency x 2 x 2)}, one quad per Rx site."""
    quads = sorted(CanonicalTargetProcessor.available_quads(path))
    out = {}
    for quad in quads:
        processor = CanonicalTargetProcessor.from_npz(path, quad=quad)
        out[quad] = processor.target_matrices()
    return out, processor.frequencies_hz


def parametric_beam_gain_dbi(beam_width_deg=60.0):
    """Directivity estimate 41253 / (theta_E theta_H) of a beam of the given -3 dB widths."""
    return 10 * np.log10(41253.0 / beam_width_deg**2)


def check_level(path, range_m, size_m, beam_width_deg):
    matrices, f = receivers(path)
    s = matrices["Rx_0/Tx"][0]
    wavelength = SPEED_OF_LIGHT / f
    sigma = 4 * np.pi * size_m**4 / (3 * wavelength**2)  # triangular trihedral, boresight
    gain = 10 ** (parametric_beam_gain_dbi(beam_width_deg) / 10)
    expected = 10 * np.log10(gain**2 * wavelength**2 * sigma / ((4 * np.pi) ** 3 * range_m**4))
    for name, (i, j) in (("HH", (0, 0)), ("VV", (1, 1))):
        level = db20(s[:, i, j])
        print(
            f"  |S_{name}|: {level.min():.2f} ... {level.max():.2f} dB over the band "
            f"(radar equation with G = {parametric_beam_gain_dbi(beam_width_deg):.1f} dBi: "
            f"{expected.mean():.2f} dB; difference {level.mean() - expected.mean():+.2f} dB)"
        )
    print(
        f"  sigma at {f[len(f) // 2] / 1e9:.1f} GHz: {10 * np.log10(sigma[len(f) // 2]):.2f} dBsm"
    )


def phase_path_m(s, f):
    """Two-way path from the phase slope: AEDT's e^{+jwt} gives arg S = -2 pi f L / c + const."""
    phase = np.unwrap(np.angle(s))
    slope = np.polyfit(f, phase, 1)[0]
    return -slope * SPEED_OF_LIGHT / (2 * np.pi)


def check_paths(path, range_m, bistatic_angle_deg):
    matrices, f = receivers(path)
    corner = np.array([range_m, 0.0, 0.0])
    tx = np.zeros(3)
    for quad, m in matrices.items():
        k = int(quad.split("/")[0].split("_")[1])
        rx = np.array([0.0, range_m * np.tan(np.radians(k * bistatic_angle_deg)), 0.0])
        expected = np.linalg.norm(rx + tx - 2 * corner)  # trihedral GO: image through corner
        naive = np.linalg.norm(corner - tx) + np.linalg.norm(corner - rx)
        for name, (i, j) in (("HH", (0, 0)), ("VV", (1, 1))):
            measured = phase_path_m(m[0, :, i, j], f)
            print(
                f"  {quad} {name}: path {measured * 1e3:9.3f} mm; GO image {expected * 1e3:9.3f}"
                f" mm (diff {(measured - expected) * 1e3:+.3f} mm); r_t + r_r "
                f"{naive * 1e3:9.3f} mm"
            )


def check_crosspol(paths):
    for path in paths:
        matrices, f = receivers(path)
        print(f"  {path.stem}:")
        for quad, m in matrices.items():
            s = m[0]
            co = np.sqrt((np.abs(s[:, 0, 0]) ** 2 + np.abs(s[:, 1, 1]) ** 2) / 2)
            for name, (i, j) in (("HV", (0, 1)), ("VH", (1, 0))):
                ratio = db20(s[:, i, j] / co)
                print(
                    f"    {quad} {name}/co: median {np.median(ratio):7.1f} dB, max "
                    f"{ratio.max():7.1f} dB; co-pol {db20(co).mean():.2f} dB"
                )


def check_strip(path):
    matrices, f = receivers(path)
    for quad, m in matrices.items():
        s = m[0]
        levels = {
            n: db20(s[:, i, j]).mean()
            for n, (i, j) in (("HH", (0, 0)), ("HV", (0, 1)), ("VH", (1, 0)), ("VV", (1, 1)))
        }
        print(f"  {quad}: " + ", ".join(f"|S_{n}| {v:.1f} dB" for n, v in levels.items()))


def check_labels(dihedral, trihedral):
    m_d, f = receivers(dihedral)
    m_t, _ = receivers(trihedral)
    s_d, s_t = m_d["Rx_0/Tx"][0], m_t["Rx_0/Tx"][0]
    for name, (i, j) in (("HH", (0, 0)), ("VV", (1, 1))):
        ratio = s_d[:, i, j] / s_t[:, i, j]
        phase = np.angle(ratio, deg=True)
        print(
            f"  S_{name}(dihedral) / S_{name}(trihedral), Rx_0: phase median "
            f"{np.median(phase):+7.1f} deg (min {phase.min():+.1f}, max {phase.max():+.1f}), "
            f"magnitude {db20(ratio).mean():+.2f} dB"
        )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("check", choices=("level", "paths", "crosspol", "strip", "labels"))
    parser.add_argument("exports", nargs="+")
    parser.add_argument("--range", type=float, default=4.0, help="Target range (m)")
    parser.add_argument("--size", type=float, default=0.04, help="Trihedral edge (m)")
    parser.add_argument("--beam-width", type=float, default=60.0, help="Beam widths (deg)")
    parser.add_argument("--bistatic-angle", type=float, default=1.0, help="bistaticAngle (deg)")
    args = parser.parse_args(argv)
    paths = [export_path(e) for e in args.exports]
    print(f"[{args.check}] {', '.join(p.name for p in paths)}")
    if args.check == "level":
        check_level(paths[0], args.range, args.size, args.beam_width)
    elif args.check == "paths":
        check_paths(paths[0], args.range, args.bistatic_angle)
    elif args.check == "crosspol":
        check_crosspol(paths)
    elif args.check == "labels":
        check_labels(paths[0], paths[1])
    else:
        check_strip(paths[0])


if __name__ == "__main__":
    main()
