"""Tests of the pure-NumPy MIMO polarimetric processing (no AEDT needed)."""

import json

import numpy as np
import pytest

from pomaar.simulator.polarimetry_processor import (
    POLARIMETRIC_CHANNELS,
    REFERENCE_DIHEDRAL,
    REFERENCE_SPHERE,
    SPEED_OF_LIGHT,
    CanonicalTargetProcessor,
    MimoPolarimetryProcessor,
    alpha_angle,
    assemble_scattering_matrix,
    calibrate_with_sphere,
    channel_powers,
    cloude_pottier,
    coherency_matrix,
    complex_from_db_deg,
    copolar_ratio,
    pauli_power_fractions,
    pauli_powers,
    pauli_vector,
    plane_wave_phase,
    polarimetric_fidelity,
    reference_dihedral,
    span,
    spatial_window,
    to_scattering_matrices,
    total_power,
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


def quad_channels(matrices):
    """(n_var, n_freq, 2, 2) -> (n_var, 4, n_freq) in the order of EXPRESSIONS."""
    return np.stack([matrices[..., r, t] for r in (0, 1) for t in (0, 1)], axis=1)


# --- per-matrix signatures ----------------------------------------------------------------------


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
    # A dihedral rolled by 45 deg is purely cross-polar and symmetric
    np.testing.assert_allclose(reference_dihedral(45.0), [[0, 1], [1, 0]], atol=1e-12)


def test_cloude_pottier_of_one_target_and_of_white_coherency():
    k = pauli_vector(np.broadcast_to(REFERENCE_SPHERE * (1 + 1j), (20, 2, 2)))
    entropy, anisotropy, mean_alpha = cloude_pottier(coherency_matrix(k))
    assert entropy == pytest.approx(0.0, abs=1e-9)
    assert mean_alpha == pytest.approx(0.0, abs=1e-6)
    entropy, anisotropy, _ = cloude_pottier(np.eye(4) / 4)
    assert entropy == pytest.approx(1.0) and anisotropy == pytest.approx(0.0)
    # Two equal mechanisms and nothing else: A = 1
    _, anisotropy, _ = cloude_pottier(np.diag([0.5, 0.5, 0.0, 0.0]))
    assert anisotropy == pytest.approx(1.0)


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


def test_sphere_calibration_restores_dihedral_phase():
    rng = np.random.default_rng(1)
    receive = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    transmit = rng.standard_normal(2) + 1j * rng.standard_normal(2)
    gains = np.outer(receive, transmit)  # g_pq = r_p t_q
    calibrated = calibrate_with_sphere(gains * REFERENCE_DIHEDRAL, gains * REFERENCE_SPHERE)
    np.testing.assert_allclose(calibrated, REFERENCE_DIHEDRAL, atol=1e-12)


def test_assemble_scattering_matrix_layout():
    s = assemble_scattering_matrix(1, 2, 3, 4)
    np.testing.assert_array_equal(s, [[1, 2], [3, 4]])


# --- (a)-(c) reconstruction and range compression -----------------------------------------------


def test_complex_from_db_deg():
    s = 0.01 * np.exp(1j * np.radians(-135.0))
    np.testing.assert_allclose(complex_from_db_deg(-40.0, -135.0), s)
    np.testing.assert_allclose(complex_from_db_deg(20 * np.log10(np.abs(s)), np.angle(s, True)), s)


def test_range_profile_peak_and_centre_frequency_phase():
    frequencies = np.linspace(77e9, 81e9, 81)
    k = 2 * np.pi * frequencies / SPEED_OF_LIGHT
    target_range = 1.2
    s = (1e-4 * np.exp(-1j * k * 2 * target_range))[None, None]
    processor = MimoPolarimetryProcessor(s, frequencies, [EXPRESSIONS[0]])
    range_m, profiles = processor.range_profiles(window="hann", oversampling=8)
    peak = np.argmax(np.abs(profiles[0, 0]))
    assert range_m[peak] == pytest.approx(target_range, abs=range_m[1])
    assert np.abs(profiles[0, 0, peak]) == pytest.approx(1e-4, rel=0.02)
    # Unambiguous range c / (2 df) = 3 m for 50 MHz steps
    assert range_m[-1] + range_m[1] == pytest.approx(SPEED_OF_LIGHT / (2 * 50e6))
    # The phase refers to the centre frequency: exp(-j k0 2R) at an exactly sampled bin
    s_exact = np.exp(-1j * k * 2 * range_m[40])[None, None]
    _, exact = MimoPolarimetryProcessor(s_exact, frequencies, [EXPRESSIONS[0]]).range_profiles(
        oversampling=8
    )
    k0 = 2 * np.pi * 79e9 / SPEED_OF_LIGHT
    assert np.angle(exact[0, 0, 40] / np.exp(-1j * k0 * 2 * range_m[40])) == pytest.approx(
        0.0, abs=1e-9
    )


# --- (d)-(g) beamforming and maps ---------------------------------------------------------------


def mimo_export(direction_deg, target_range=1.5, scattering=REFERENCE_SPHERE):
    """
    A synthetic 2D MIMO array observing one point target with scattering matrix `scattering`
    at (azimuth, elevation) `direction_deg`: 4 Rx H/V pairs at x = 0..3 d and 4 Tx H/V pairs at
    (0 | 4 d, 0 | d), d = lambda0 / 2, i.e. a filled 8 x 2 virtual array of 16 channels per
    polarimetric channel (no grating lobes).
    """
    frequencies = np.linspace(77e9, 81e9, 81)
    k = 2 * np.pi * frequencies / SPEED_OF_LIGHT
    az, el = np.radians(direction_deg)
    r = np.array([np.cos(el) * np.sin(az), np.sin(el), np.cos(el) * np.cos(az)])
    d = SPEED_OF_LIGHT / 79e9 / 2
    tx = {i: np.array([4 * d * (i % 2), d * (i // 2), 0.0]) for i in range(4)}
    rx = {j: np.array([j * d, 0.0, 0.0]) for j in range(4)}
    expressions, channels, positions = [], [], {}
    for t, p_tx in tx.items():
        for pol_t in "HV":
            positions[f"Tx_{t}_{pol_t}"] = p_tx
            for j, p_rx in rx.items():
                for pol_r in "HV":
                    positions[f"Rx_{j}_{pol_r}"] = p_rx
                    expressions.append(f"S(A1_Port_Rx_{j}_{pol_r},A1_Port_Tx_{t}_{pol_t})")
                    value = scattering["HV".index(pol_r), "HV".index(pol_t)]
                    path = 2 * target_range - (p_tx + p_rx) @ r
                    channels.append(value * np.exp(-1j * k * path))
    s = np.array(channels)[None]
    return s, frequencies, expressions, positions


def test_beamforming_peaks_at_the_target_direction():
    s, frequencies, expressions, positions = mimo_export((12.0, -6.0))
    processor = MimoPolarimetryProcessor(s, frequencies, expressions, positions=positions)
    assert processor.positions_known
    assert {q: len(i) for q, i in processor.polarimetric_channels.items()} == {
        q: 16 for q in POLARIMETRIC_CHANNELS
    }
    azimuth, elevation = np.arange(-30, 31, 1.0), np.arange(-20, 21, 1.0)
    maps = processor.azimuth_elevation_maps(azimuth, elevation)
    power = total_power(maps["channels"])[0]
    peak_az, peak_el = np.unravel_index(np.argmax(power), power.shape)
    assert (azimuth[peak_az], elevation[peak_el]) == (12.0, -6.0)
    # A unique main lobe: away from it, every direction is at least 10 dB down
    far = np.abs(azimuth[:, None] - 12.0) > 15
    assert np.max(power[far[:, 0]]) < 0.1 * power.max()
    assert maps["range_m"][0] == pytest.approx(1.5, abs=0.02)
    # Range-azimuth map at the target elevation peaks at the same azimuth and range
    ra = processor.range_azimuth_maps(azimuth, elevation_deg=-6.0)
    ra_power = total_power(ra["channels"])[0]
    r_bin, a_bin = np.unravel_index(np.argmax(ra_power), ra_power.shape)
    assert azimuth[a_bin] == 12.0 and ra["range_m"][r_bin] == pytest.approx(1.5, abs=0.02)


def test_map_quantities_of_a_dihedral():
    s, frequencies, expressions, positions = mimo_export((0.0, 0.0), scattering=REFERENCE_DIHEDRAL)
    processor = MimoPolarimetryProcessor(s, frequencies, expressions, positions=positions)
    maps = processor.azimuth_elevation_maps([0.0], [0.0])
    powers = channel_powers(maps["channels"])
    assert powers["HH"][0, 0, 0] == pytest.approx(2.0)  # relative to the total power (mean of 4)
    assert powers["HV"][0, 0, 0] == pytest.approx(0.0, abs=1e-20)
    pauli = pauli_powers(maps["channels"])[0, 0, 0]
    assert pauli[1] / pauli.sum() == pytest.approx(1.0)  # all power in HH - VV


def test_spatial_window_tapers_the_aperture_edges():
    positions = np.array([[[x, 0.0, 0.0] for x in range(5)]], dtype=float)
    weights = spatial_window(positions, "hann")[0]
    assert np.all(weights > 0) and weights[2] == pytest.approx(1.0)
    np.testing.assert_allclose(weights, weights[::-1])
    assert weights[0] < weights[1] < weights[2]
    np.testing.assert_allclose(spatial_window(positions, None), 1.0)


def test_raw_forwards_the_export_untouched():
    matrices = random_matrices((2, len(FREQUENCIES)))
    processor = MimoPolarimetryProcessor(quad_channels(matrices), FREQUENCIES, EXPRESSIONS)
    processor.inject_noise(10.0, rng=0)
    np.testing.assert_array_equal(processor.raw()["s"], quad_channels(matrices))
    assert processor.raw()["expressions"] == EXPRESSIONS


# --- canonical targets --------------------------------------------------------------------------


def test_target_matrices_at_boresight_are_the_exported_quad():
    matrices = random_matrices((2, len(FREQUENCIES)))
    processor = CanonicalTargetProcessor(quad_channels(matrices), FREQUENCIES, EXPRESSIONS)
    assert not processor.positions_known  # no port table: positions 0, exact at boresight
    np.testing.assert_allclose(processor.target_matrices(), matrices)


def test_target_matrices_compensate_the_linear_term_off_boresight():
    s, frequencies, expressions, positions = mimo_export(
        (8.0, 3.0), scattering=reference_dihedral(30)
    )
    processor = CanonicalTargetProcessor(s, frequencies, expressions, positions=positions)
    matrices = processor.target_matrices(8.0, 3.0)
    expected = (
        reference_dihedral(30)
        * np.exp(-1j * 2 * np.pi * frequencies / SPEED_OF_LIGHT * 3.0)[:, None, None]
    )
    np.testing.assert_allclose(matrices[0], expected, atol=1e-12)
    assert set(processor.bounce_parity(matrices).ravel()) == {"even"}
    assert np.min(processor.fidelity(matrices, reference_dihedral(30))) == pytest.approx(1.0)
    features = processor.signatures(matrices=matrices)
    np.testing.assert_allclose(features["copol_phase_difference_deg"] % 360, 180.0)
    np.testing.assert_allclose(features["entropy"], 0.0, atol=1e-9)


def test_canonical_features_tell_odd_and_even_bounces_apart():
    sphere = 1e-3 * np.exp(1j * 0.4) * REFERENCE_SPHERE
    leaky = sphere + np.array([[0, 1e-6], [3e-6, 0]])
    processor = CanonicalTargetProcessor
    assert processor.bounce_parity(sphere) == "odd"
    assert processor.bounce_parity(REFERENCE_DIHEDRAL) == "even"
    co_h, co_v = processor.copol_crosspol_difference_db(leaky)
    assert co_h == pytest.approx(20 * np.log10(1e-3 / 3e-6))
    assert co_v == pytest.approx(60.0)


def test_signatures_are_target_agnostic_quantities(tmp_path):
    sphere = np.zeros((2, 4, len(FREQUENCIES)), dtype=complex)
    sphere[:, [0, 3]] = 1e-5 * np.exp(-1j * WAVENUMBER * 5.85)
    sphere[:, 2] = 1e-8  # a little VH
    np.savez(
        tmp_path / "sphere.npz",
        s=sphere,
        frequencies_hz=FREQUENCIES,
        sweep_values=np.array([1.41, 2.41]),
        expressions=np.asarray(EXPRESSIONS),
        metadata=json.dumps({"setup_sweep": "Setup : Sweep_Scattered"}),
    )
    processor = CanonicalTargetProcessor.from_npz(tmp_path / "sphere.npz")
    features = processor.signatures()
    assert not any("fidelity" in key or "reference" in key for key in features)
    np.testing.assert_allclose(features["power_hh_db"], -100.0)
    np.testing.assert_allclose(features["power_vh_db"], -160.0)
    np.testing.assert_allclose(features["copol_ratio_db"], 0.0, atol=1e-9)
    np.testing.assert_allclose(features["co_minus_cross_txh_db"], 60.0)
    assert features["entropy"].shape == (2,) and features["pauli_fractions_db"].shape[-1] == 4


# --- builder exports ----------------------------------------------------------------------------


def builder_export(tmp_path, matrices, ports, rows=None):
    """Writes an export whose expressions use builder port names (Rx first, then Tx)."""
    expressions, channels = [], []
    for rx_index, rx_pol in enumerate("HV"):
        for tx_index, tx_pol in enumerate("HV"):
            expressions.append(f"S({ports['Rx' + rx_pol]},{ports['Tx' + tx_pol]})")
            channels.append(matrices[..., rx_index, tx_index])
    metadata = {"setup_sweep": "Setup : Sweep_Scattered"}
    if rows is not None:
        metadata["ports"] = rows
    path = tmp_path / "builder.npz"
    np.savez(
        path,
        s=np.stack(channels, axis=1)[:, :, :],
        frequencies_hz=FREQUENCIES[: matrices.shape[1]],
        sweep_values=np.arange(matrices.shape[0], dtype=float),
        expressions=np.asarray(expressions),
        metadata=json.dumps(metadata),
    )
    return path


CLUSTER_PORTS = {
    "RxH": "Cluster1_Port_Rx_H",
    "RxV": "Cluster1_Port_Rx_V",
    "TxH": "Cluster1_Port_Tx_H",
    "TxV": "Cluster1_Port_Tx_V",
}


def test_builder_port_names_map_like_hand_built_ones():
    s = np.arange(4)[None, :, None] * np.ones((1, 4, 3))
    expressions = [
        f"S({CLUSTER_PORTS['Rx' + r]},{CLUSTER_PORTS['Tx' + t]})" for r in "HV" for t in "HV"
    ]
    np.testing.assert_array_equal(to_scattering_matrices(s, expressions)[0, 0], [[0, 1], [2, 3]])
    # Expression order does not matter
    np.testing.assert_array_equal(
        to_scattering_matrices(s[:, ::-1], expressions[::-1])[0, 0], [[0, 1], [2, 3]]
    )


def test_quad_selection_keeps_one_quad(tmp_path):
    matrices = random_matrices((1, 3))
    ports = {key: f"A1_Port_{key[:2]}_1_{key[2]}" for key in CLUSTER_PORTS}
    other = {key: f"A1_Port_{key[:2]}_2_{key[2]}" for key in CLUSTER_PORTS}
    path = builder_export(tmp_path, matrices, ports)
    # Add the channels of a second Tx element so that two quads exist
    data = dict(np.load(path))
    extra = [f"S({ports['Rx' + r]},{other['Tx' + t]})" for r in "HV" for t in "HV"]
    data["expressions"] = np.concatenate([data["expressions"], extra])
    data["s"] = np.concatenate([data["s"], 2 * data["s"]], axis=1)
    np.savez(path, **data)
    assert MimoPolarimetryProcessor.available_quads(path) == ["Rx_1/Tx_1", "Rx_1/Tx_2"]
    everything = MimoPolarimetryProcessor.from_npz(path)
    assert {q: len(i) for q, i in everything.polarimetric_channels.items()} == {
        q: 2 for q in POLARIMETRIC_CHANNELS
    }
    processor = CanonicalTargetProcessor.from_npz(path, quad="Rx_1/Tx_2")
    np.testing.assert_allclose(processor.target_matrices(), 2 * matrices)


def test_port_table_signs_and_positions(tmp_path):
    # pi_Tx = pi_Rx = -1 (the old polarimetric_ula yaws): HV and VH are measured negated
    signs = {"RxH": -1, "RxV": 1, "TxH": 1, "TxV": -1}
    truth = random_matrices((2, 3))
    measured = truth * np.outer([signs["RxH"], signs["RxV"]], [signs["TxH"], signs["TxV"]])
    rows = [
        {
            "excitation": CLUSTER_PORTS[key],
            "port": CLUSTER_PORTS[key].removeprefix("Cluster1_"),
            "sign": signs[key],
            "phase_centre_mm": [float(i), 2.0, 0.0],
        }
        for i, key in enumerate(CLUSTER_PORTS)
    ]
    path = builder_export(tmp_path, measured, CLUSTER_PORTS, rows)
    processor = CanonicalTargetProcessor.from_npz(path)
    assert processor.signs_corrected and processor.positions_known
    np.testing.assert_allclose(processor.target_matrices(), truth)
    raw = CanonicalTargetProcessor.from_npz(path, correct_signs=False)
    np.testing.assert_allclose(raw.target_matrices(), measured)
    # Virtual position of HV = TxV (3 mm, 2 mm) + RxH (0, 2 mm)
    np.testing.assert_allclose(processor.virtual_positions("HV")[0, 0], [0.003, 0.004, 0.0])
    # Without a port table there is nothing to correct or place
    bare = MimoPolarimetryProcessor.from_npz(builder_export(tmp_path, measured, CLUSTER_PORTS))
    assert not bare.signs_corrected and not bare.positions_known
