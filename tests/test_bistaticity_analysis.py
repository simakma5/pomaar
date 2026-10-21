"""Tests of the per-receiver pooling of `bistaticity_analysis` (no AEDT needed)."""

import json

import numpy as np
import pytest

from pomaar.simulator.bistaticity_analysis import pool_receivers, steering_direction

FREQUENCIES = np.linspace(74e9, 84e9, 11)
SITES = ("Rx_0", "Rx_1", "Rx_2")


def receiver_export(path, angles_deg, range_m=4.0):
    """A sphere-like export of one Tx and three Rx sites at k * angle, with scene geometry."""
    expressions, channels = [], []
    for site in SITES:
        for rx in "HV":
            for tx in "HV":
                expressions.append(f"S({site}_{rx}1_p1,Tx_{tx}1_p1)")
                channels.append(1e-5 if rx == tx else 1e-9)
    s = np.tile(np.asarray(channels)[None, :, None], (len(angles_deg), 1, len(FREQUENCIES)))
    positions = {"Tx_H1_p1": [], "Tx_V1_p1": []}
    for angle in angles_deg:
        positions["Tx_H1_p1"].append([0.0, 0.0, 0.0])
        positions["Tx_V1_p1"].append([0.0, 0.0, 0.0])
        for k, site in enumerate(SITES):
            for pol in "HV":
                y = range_m * np.tan(np.radians(k * angle))
                positions.setdefault(f"{site}_{pol}1_p1", []).append([0.0, y, 0.0])
    metadata = {
        "setup_sweep": "Setup : Sweep_Scattered",
        "sweep_variable": "bistaticAngle",
        "positions_m": positions,
        "target_m": [[range_m, 0.0, 0.0]] * len(angles_deg),
    }
    np.savez(
        path,
        s=s.astype(complex),
        frequencies_hz=FREQUENCIES,
        sweep_values=np.asarray(angles_deg, dtype=float),
        expressions=np.asarray(expressions),
        metadata=json.dumps(metadata),
    )
    return path


def test_pool_receivers_sorts_by_bistatic_angle_and_keeps_equal_angles_once(tmp_path):
    path = receiver_export(tmp_path / "pooled.npz", [1.0, 2.0])
    features, angles, receivers, values = pool_receivers(path)
    # Rx_0 is monostatic at both variations, and Rx_1 at 2 deg sits where Rx_2 at 1 deg does
    np.testing.assert_allclose(angles, [0.0, 1.0, 2.0, 4.0], atol=1e-6)
    assert receivers == ["Rx_0", "Rx_1", "Rx_1", "Rx_2"]
    assert features["ldr_txh_db"].shape == (4, len(FREQUENCIES))
    np.testing.assert_allclose(features["ldr_txh_db"], -80.0)
    assert values[1] == pytest.approx(1.0)


def test_steering_direction_defaults_to_a_fixed_target_direction():
    fixed = {"target_steering_deg": [[10.0, 5.0], [10.0, 5.0]]}
    assert steering_direction(fixed)[:2] == (10.0, 5.0)
    assert steering_direction(fixed, azimuth=0.0)[:2] == (0.0, 5.0)
    assert steering_direction(fixed, 1.0, 2.0)[:2] == (1.0, 2.0)


def test_steering_direction_falls_back_to_boresight():
    assert steering_direction({})[:2] == (0.0, 0.0)
    moving = {"target_steering_deg": [[0.0, 0.0], [10.0, 0.0]]}
    azimuth, elevation, source = steering_direction(moving)
    assert (azimuth, elevation) == (0.0, 0.0) and "moves" in source
