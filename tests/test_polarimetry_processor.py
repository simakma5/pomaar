"""Tests of the pure-NumPy polarimetric routines (no AEDT needed)."""

import json

import numpy as np
import pytest

from pomaar.simulator.polarimetry_processor import (
    REFERENCE_DIHEDRAL,
    REFERENCE_SPHERE,
    SPEED_OF_LIGHT,
    MimoPolarimetryProcessor,
    alpha_angle,
    assemble_scattering_matrix,
    coherency_matrix,
    copolar_ratio,
    entropy_alpha,
    gate_valid_band,
    pauli_power_fractions,
    pauli_vector,
    plane_wave_phase,
    polarimetric_fidelity,
    span,
    spherical_wave_phase,
    time_gate,
    to_scattering_matrices,
)

EXPRESSIONS = [
    "S(HornCluster1_RxH,HornCluster1_TxH)",
    "S(HornCluster1_RxH,HornCluster1_TxV)",
    "S(HornCluster1_RxV,HornCluster1_TxH)",
    "S(HornCluster1_RxV,HornCluster1_TxV)",
]
FREQUENCIES = np.linspace(74e9, 84e9, 201)
WAVENUMBER = 2 * np.pi * FREQUENCIES / SPEED_OF_LIGHT


def random_matrices(shape, seed=0):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((*shape, 2, 2)) + 1j * rng.standard_normal((*shape, 2, 2))


def cluster(spacing_m):
    """The DualPolHornCluster square: co-polar pairs on the diagonals, spacing apart."""
    offset = spacing_m / (2 * np.sqrt(2))
    return {
        "TxV": np.array([offset, offset, 0.0]),
        "RxV": np.array([-offset, -offset, 0.0]),
        "TxH": np.array([offset, -offset, 0.0]),
        "RxH": np.array([-offset, offset, 0.0]),
    }


def test_aedt_expressions_map_to_receive_transmit_indices():
    s = np.arange(4)[None, :, None] * np.ones((1, 4, 3))
    matrices = to_scattering_matrices(s, EXPRESSIONS)
    # S(RxH,TxV) is S_HV = S[rx=H, tx=V]
    np.testing.assert_array_equal(matrices[0, 0], [[0, 1], [2, 3]])
    with pytest.raises(ValueError):
        to_scattering_matrices(s, EXPRESSIONS[:3] + EXPRESSIONS[:1])


def test_pauli_vector_preserves_span_and_fractions_sum_to_one():
    s = random_matrices((5, 7))
    np.testing.assert_allclose(np.sum(np.abs(pauli_vector(s)) ** 2, axis=-1), span(s))
    np.testing.assert_allclose(pauli_power_fractions(s).sum(axis=-1), 1.0)


def test_antisymmetric_component_vanishes_for_reciprocal_matrices():
    s = random_matrices((10,))
    s[..., 1, 0] = s[..., 0, 1]
    np.testing.assert_allclose(pauli_vector(s)[..., 3], 0.0, atol=1e-12)


def test_fidelity_and_alpha_of_canonical_targets():
    scale = 3.0 * np.exp(1j * 0.7)
    assert polarimetric_fidelity(scale * REFERENCE_SPHERE, REFERENCE_SPHERE) == pytest.approx(1.0)
    assert polarimetric_fidelity(REFERENCE_DIHEDRAL, REFERENCE_SPHERE) == pytest.approx(0.0)
    assert alpha_angle(scale * REFERENCE_SPHERE) == pytest.approx(0.0, abs=1e-6)
    assert alpha_angle(REFERENCE_DIHEDRAL) == pytest.approx(90.0)
    assert np.angle(copolar_ratio(REFERENCE_DIHEDRAL), deg=True) == pytest.approx(180.0)


def test_entropy_is_zero_for_one_target_and_one_for_white_coherency():
    k = pauli_vector(np.broadcast_to(REFERENCE_SPHERE * (1 + 1j), (20, 2, 2)))
    entropy, mean_alpha = entropy_alpha(coherency_matrix(k))
    assert entropy == pytest.approx(0.0, abs=1e-9)
    assert mean_alpha == pytest.approx(0.0, abs=1e-6)
    entropy, _ = entropy_alpha(np.eye(4) / 4)
    assert entropy == pytest.approx(1.0)


@pytest.mark.parametrize("coupling_db", [-35.6, -67.4])
def test_time_gate_recovers_target_next_to_strong_coupling(coupling_db):
    # Measured situation: sphere at 5.851 m (wraps to -0.145 m), coupling at +0.056 m
    sphere = 10 ** (-97.5 / 20) * np.exp(-1j * WAVENUMBER * 5.851)
    coupling = 10 ** (coupling_db / 20) * np.exp(-1j * WAVENUMBER * 0.056)
    window = ("kaiser", 14)
    gated = time_gate(sphere + coupling, FREQUENCIES, 5.851, 0.06, window=window)
    band = gate_valid_band(FREQUENCIES, window=window)
    error = gated[band] / sphere[band]
    assert np.max(np.abs(20 * np.log10(np.abs(error)))) < 0.05
    assert np.max(np.abs(np.angle(error, deg=True))) < 0.3


def test_linear_term_between_cross_polar_channels():
    spacing = 8.41 * 3.795e-3
    positions = cluster(spacing)
    theta = np.radians(2.0)
    direction = np.array([0.0, np.sin(theta), np.cos(theta)])
    hv = plane_wave_phase(FREQUENCIES, positions["TxV"], positions["RxH"], direction)
    vh = plane_wave_phase(FREQUENCIES, positions["TxH"], positions["RxV"], direction)
    expected = np.sqrt(2) * WAVENUMBER * spacing * np.sin(theta)
    np.testing.assert_allclose(np.unwrap(np.angle(hv / vh)), expected, rtol=1e-9)
    # The co-polar channels share their virtual element at the centroid
    hh = plane_wave_phase(FREQUENCIES, positions["TxH"], positions["RxH"], direction)
    np.testing.assert_allclose(hh, 1.0)


def test_spherical_compensation_removes_all_geometric_phase():
    positions = cluster(5.41 * 3.795e-3)
    target = np.array([0.0, 0.3, 0.8])  # near and off boresight
    matrices = np.zeros((1, len(FREQUENCIES), 2, 2), dtype=complex)
    for rx, rx_name in enumerate("HV"):
        for tx, tx_name in enumerate("HV"):
            matrices[0, :, rx, tx] = spherical_wave_phase(
                FREQUENCIES, positions[f"Tx{tx_name}"], positions[f"Rx{rx_name}"], target
            )
    processor = MimoPolarimetryProcessor(matrices, FREQUENCIES)
    processor.compensate_geometry(positions, target)
    np.testing.assert_allclose(processor.scattering_matrices, 1.0, atol=1e-9)


def test_sphere_calibration_restores_dihedral_phase():
    rng = np.random.default_rng(1)
    receive = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    transmit = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    gains = np.outer(receive, transmit)  # g_pq = r_p t_q
    sphere = MimoPolarimetryProcessor((gains * REFERENCE_SPHERE)[None, None], FREQUENCIES[:1])
    dihedral = MimoPolarimetryProcessor((gains * REFERENCE_DIHEDRAL)[None, None], FREQUENCIES[:1])
    dihedral.calibrate_with_sphere(sphere)
    np.testing.assert_allclose(dihedral.scattering_matrices[0, 0], REFERENCE_DIHEDRAL, atol=1e-12)


def test_from_npz_subtracts_background_and_signatures_of_ideal_sphere(tmp_path):
    sphere = np.zeros((2, 4, len(FREQUENCIES)), dtype=complex)
    sphere[:, [0, 3]] = 1e-5 * np.exp(-1j * WAVENUMBER * 5.85)
    coupling = 1e-2 * np.exp(-1j * WAVENUMBER * 0.06) * np.ones((2, 4, 1))
    common = {
        "frequencies_hz": FREQUENCIES,
        "sweep_values": np.array([1.41, 2.41]),
        "expressions": np.asarray(EXPRESSIONS),
        "metadata": json.dumps({"setup_sweep": "Setup : Sweep"}),
    }
    np.savez(tmp_path / "total.npz", s=sphere + coupling, **common)
    np.savez(tmp_path / "incident.npz", s=coupling, **common)

    processor = MimoPolarimetryProcessor.from_npz(tmp_path / "total.npz", tmp_path / "incident.npz")
    np.testing.assert_allclose(
        processor.scattering_matrices, to_scattering_matrices(sphere, EXPRESSIONS), atol=1e-15
    )
    signatures = processor.signatures(REFERENCE_SPHERE)
    np.testing.assert_allclose(signatures["copol_ratio"], 1.0)
    np.testing.assert_allclose(signatures["infidelity"], 0.0, atol=1e-12)
    np.testing.assert_allclose(signatures["entropy_bistatic"], 0.0, atol=1e-9)


def test_assemble_scattering_matrix_layout():
    s = assemble_scattering_matrix(1, 2, 3, 4)
    np.testing.assert_array_equal(s, [[1, 2], [3, 4]])
