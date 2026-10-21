#!/usr/bin/env python3
"""
S3/S4 (host side): steering consistency between the SBR+ scene and the processor. A dihedral at
roll 45 deg sits 10 deg in azimuth (towards H) and 5 deg in elevation (towards V) off boresight,
in an x-boresight scene and in a z-boresight scene. The cluster's HV and VH virtual channels lie
on opposite sides of the H/V diagonal, so the cross-pol phase arg(S_HV S_VH*) of a reciprocal
target is ~0 only when the channels are steered at the target. Printed per export and steering
direction (centre frequency and band maximum).

    python -m pomaar.simulator.investigation.validation.steering_check X.npz:az,el ... Z.npz:...
"""

import sys

import numpy as np

from pomaar.simulator.polarimetry_processor import CanonicalTargetProcessor


def main(argv=None):
    for item in argv or sys.argv[1:]:
        path, *directions = item.split(":")
        processor = CanonicalTargetProcessor.from_npz(path)
        print(
            f"\n{path} (positions known: {processor.positions_known}, "
            f"signs corrected: {processor.signs_corrected})"
        )
        centre = len(processor.frequencies_hz) // 2
        for direction in directions:
            az, el = (float(v) for v in direction.split(","))
            s = processor.target_matrices(az, el)[0]  # first variation, (n_freq, 2, 2)
            phase = np.degrees(np.angle(s[:, 0, 1] * np.conj(s[:, 1, 0])))
            ratio = 20 * np.log10(np.abs(s[:, 0, 1] / s[:, 1, 0]))
            print(
                f"  steered to ({az:g}, {el:g}): arg(S_HV S_VH*) {phase[centre]:+7.2f} deg at "
                f"f0 (band max |.| {np.abs(phase).max():6.2f}), |S_HV/S_VH|^2 "
                f"{ratio[centre]:+.3f} dB"
            )


if __name__ == "__main__":
    main()
