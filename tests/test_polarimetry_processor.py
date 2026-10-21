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
    IncompleteSweepError,
    MimoPolarimetryProcessor,
    alpha_angle,
    assemble_scattering_matrix,
    bistatic_angle_deg,
    calibrate_with_sphere,
    channel_powers,
    check_complete_sweep,
    cloude_pottier,
    coherency_matrix,
    complex_from_db_deg,
    copolar_ratio,
    identical_variations,
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
from pomaar.simulator.scene_frames import array_to_scene

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


def mimo_export(direction_deg, target_range=1.5, scattering=REFERENCE_SPHERE, boresight="x"):
    """
    A synthetic 2D MIMO array observing one point target with scattering matrix `scattering`
    at (azimuth, elevation) `direction_deg`: 4 Rx H/V pairs at H = 0..3 d and 4 Tx H/V pairs at
    (H, V) = (0 | 4 d, 0 | d), d = lambda0 / 2, i.e. a filled 8 x 2 virtual array of 16 channels
    per polarimetric channel (no grating lobes). Built in the array frame (H = x, V = y); the
    positions are returned in the scene frame of `boresight`.
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
    positions = {port: array_to_scene(p, boresight) for port, p in positions.items()}
    return s, frequencies, expressions, positions


@pytest.mark.parametrize("boresight", ["x", "z"])
def test_beamforming_peaks_at_the_target_direction(boresight):
    s, frequencies, expressions, positions = mimo_export((12.0, -6.0), boresight=boresight)
    processor = MimoPolarimetryProcessor(
        s, frequencies, expressions, positions=positions, boresight=boresight
    )
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


def test_beamforming_is_the_same_in_both_frames():
    # The same export, with positions in either frame, gives the same maps and spatial weights
    maps = {}
    for boresight in ("x", "z"):
        s, frequencies, expressions, positions = mimo_export((7.0, 4.0), boresight=boresight)
        processor = MimoPolarimetryProcessor(
            s, frequencies, expressions, positions=positions, boresight=boresight
        )
        maps[boresight] = processor.azimuth_elevation_maps(
            np.arange(-20, 21, 2.0), np.arange(-10, 11, 2.0), spatial_window="hann"
        )
        assert processor.spatial_weights("HH", "hann").std() > 0  # a taper over both aperture axes
    for q in POLARIMETRIC_CHANNELS:
        np.testing.assert_allclose(maps["x"]["channels"][q], maps["z"]["channels"][q])


def test_map_quantities_of_a_dihedral():
    s, frequencies, expressions, positions = mimo_export((0.0, 0.0), scattering=REFERENCE_DIHEDRAL)
    processor = MimoPolarimetryProcessor(s, frequencies, expressions, positions=positions)
    maps = processor.azimuth_elevation_maps([0.0], [0.0])
    powers = channel_powers(maps["channels"])
    assert powers["HH"][0, 0, 0] == pytest.approx(2.0)  # relative to the total power (mean of 4)
    assert powers["HV"][0, 0, 0] == pytest.approx(0.0, abs=1e-20)
    pauli = pauli_powers(maps["channels"])[0, 0, 0]
    assert pauli[1] / pauli.sum() == pytest.approx(1.0)  # all power in HH - VV


@pytest.mark.parametrize("boresight", ["x", "z"])
def test_spatial_window_tapers_the_aperture_edges(boresight):
    # Five channels along H (array frame x), in the scene frame
    positions = array_to_scene([[[h, 0.0, 0.0] for h in range(5)]], boresight)
    weights = spatial_window(positions, "hann", boresight)[0]
    assert np.all(weights > 0) and weights[2] == pytest.approx(1.0)
    np.testing.assert_allclose(weights, weights[::-1])
    assert weights[0] < weights[1] < weights[2]
    np.testing.assert_allclose(spatial_window(positions, None, boresight), 1.0)
    # Spread along boresight instead, the channels are not tapered
    along_boresight = array_to_scene([[[0.0, 0.0, z] for z in range(5)]], boresight)
    np.testing.assert_allclose(spatial_window(along_boresight, "hann", boresight), 1.0)


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


@pytest.mark.parametrize("boresight", ["x", "z"])
def test_target_matrices_compensate_the_linear_term_off_boresight(boresight):
    s, frequencies, expressions, positions = mimo_export(
        (8.0, 3.0), scattering=reference_dihedral(30), boresight=boresight
    )
    processor = CanonicalTargetProcessor(
        s, frequencies, expressions, positions=positions, boresight=boresight
    )
    matrices = processor.target_matrices(8.0, 3.0)
    expected = (
        reference_dihedral(30)
        * np.exp(-1j * 2 * np.pi * frequencies / SPEED_OF_LIGHT * 3.0)[:, None, None]
    )
    np.testing.assert_allclose(matrices[0], expected, atol=1e-12)
    assert set(processor.bounce_parity(matrices).ravel()) == {"even"}
    assert np.min(processor.fidelity(matrices, reference_dihedral(30))) == pytest.approx(1.0)
    features = processor.signatures(matrices=matrices)
    np.testing.assert_allclose(features["hh_vv_phase_deg"] % 360, 180.0)
    np.testing.assert_allclose(features["entropy"], 0.0, atol=1e-9)


def test_canonical_features_tell_odd_and_even_bounces_apart():
    sphere = 1e-3 * np.exp(1j * 0.4) * REFERENCE_SPHERE
    leaky = sphere + np.array([[0, 1e-6], [3e-6, 0]])
    processor = CanonicalTargetProcessor
    assert processor.bounce_parity(sphere) == "odd"
    assert processor.bounce_parity(REFERENCE_DIHEDRAL) == "even"
    ldr_h, ldr_v = processor.ldr_db(leaky)
    assert ldr_h == pytest.approx(20 * np.log10(3e-6 / 1e-3))
    assert ldr_v == pytest.approx(-60.0)


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
    np.testing.assert_allclose(features["zdr_db"], 0.0, atol=1e-9)
    np.testing.assert_allclose(features["ldr_txh_db"], -60.0)
    assert features["entropy"].shape == (2,) and features["pauli_fractions_db"].shape[-1] == 4


def test_co_cross_signatures_pair_the_receivers_of_one_transmission(tmp_path):
    # S_HH = 1, S_HV = 0.1j, S_VH = 0.01, S_VV = -1, times a common path phase
    channels = np.array([1.0, 0.1j, 0.01, -1.0])
    s = 1e-5 * channels[None, :, None] * np.exp(-1j * WAVENUMBER * 5.85)[None, None, :]
    np.savez(
        tmp_path / "leaky.npz",
        s=s,
        frequencies_hz=FREQUENCIES,
        sweep_values=np.array([1.0]),
        expressions=np.asarray(EXPRESSIONS),
        metadata=json.dumps({"setup_sweep": "Setup : Sweep_Scattered"}),
    )
    features = CanonicalTargetProcessor.from_npz(tmp_path / "leaky.npz").signatures()
    # LDR = |S_cross|^2 / |S_co|^2 and arg(S_cross S_co*): Tx H (S_VH / S_HH = 0.01), Tx V
    # (S_HV / S_VV = -0.1j)
    expected = {"txh": (-40.0, 0.0), "txv": (-20.0, -90.0)}
    for side, (ratio_db, phase_deg) in expected.items():
        np.testing.assert_allclose(features[f"ldr_{side}_db"], ratio_db, atol=1e-9)
        phase = features[f"cross_co_phase_{side}_deg"]
        np.testing.assert_allclose(np.cos(np.radians(phase - phase_deg)), 1.0, atol=1e-9)
    # Z_DR = |S_HH|^2 / |S_VV|^2 and arg(S_HH S_VV*) of S_HH = 1, S_VV = -1
    np.testing.assert_allclose(features["zdr_db"], 0.0, atol=1e-9)
    np.testing.assert_allclose(np.abs(features["hh_vv_phase_deg"]), 180.0)


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
    # Virtual position of HV = TxV (3 mm, 2 mm) + RxH (0, 2 mm) in the array frame (H, V,
    # boresight), i.e. (boresight, H, V) in the default ISO frame
    np.testing.assert_allclose(processor.virtual_positions("HV")[0, 0], [0.0, 0.003, 0.004])
    in_array_frame = CanonicalTargetProcessor.from_npz(path, boresight="z")
    np.testing.assert_allclose(in_array_frame.virtual_positions("HV")[0, 0], [0.003, 0.004, 0.0])
    # Without a port table there is nothing to correct or place
    bare = MimoPolarimetryProcessor.from_npz(builder_export(tmp_path, measured, CLUSTER_PORTS))
    assert not bare.signs_corrected and not bare.positions_known


# --- incomplete sweeps --------------------------------------------------------------------------


def test_incomplete_sweeps_are_refused_with_the_variations_named():
    frequencies = np.linspace(74e9, 84e9, 1001)
    s = np.ones((3, 4, len(frequencies)), dtype=complex)
    # Variations 0 and 2 were solved on the old 201-point grid: NaN elsewhere (as AEDT returns)
    old_grid = np.isin(np.arange(len(frequencies)), np.arange(0, 1001, 5))
    s[[0, 2]] = np.where(old_grid, 1.0, np.nan)
    with pytest.raises(IncompleteSweepError) as error:
        check_complete_sweep(s, frequencies, [1.41, 1.66, 2.41], "copolarSpacingLambda")
    message = str(error.value)
    assert "2 of 3 values of copolarSpacingLambda" in message
    assert "copolarSpacingLambda = 1.41: 201 of 1001 frequencies (step 50 MHz)" in message
    assert "complete: 1.66" in message
    with pytest.raises(IncompleteSweepError):
        MimoPolarimetryProcessor(s, frequencies, EXPRESSIONS, [1.41, 1.66, 2.41])
    check_complete_sweep(s[[1]], frequencies, [1.66])  # a complete export passes


def test_unsolved_variation_and_inf_are_incomplete():
    s = np.ones((2, 4, 5), dtype=complex)
    s[1] = np.nan
    s[0, 2, 3] = np.inf
    with pytest.raises(IncompleteSweepError, match="no data at all"):
        check_complete_sweep(s, np.linspace(1e9, 2e9, 5), None)


def test_identical_variations_are_named():
    rng = np.random.default_rng(1)
    s = rng.standard_normal((3, 4, 5)) + 1j * rng.standard_normal((3, 4, 5))
    assert identical_variations(s, [1.0, 2.0, 3.0], "spacing") is None
    same = np.repeat(s[:1], 3, axis=0)
    assert "All 3 values of spacing give identical data" in identical_variations(
        same, [1.0, 2.0, 3.0], "spacing"
    )
    s[2] = s[0]
    assert "(1, 3)" in identical_variations(s, [1.0, 2.0, 3.0], "spacing")
    assert identical_variations(s[:1], [1.0]) is None


# --- analysis plots -----------------------------------------------------------------------------


def test_phase_bands_are_taken_around_the_line():
    from pomaar.simulator.bistaticity_analysis import band_limits

    # A phase scattered around +-180 deg: plain min-max spans the circle, the wrapped one does not
    phase = np.array([[178.0, -179.0, 179.5, -177.0]])
    low, high = band_limits(phase, 2, circular=True)
    np.testing.assert_allclose([low[0], high[0]], [178.0, 183.0])
    low, high = band_limits(phase, 2)
    np.testing.assert_allclose([low[0], high[0]], [-179.0, 179.5])


def test_analysis_writes_its_figures_and_warns_on_identical_variations(tmp_path, capsys):
    from pomaar.simulator.bistaticity_analysis import main

    dihedral = np.zeros((3, 4, len(FREQUENCIES)), dtype=complex)
    dihedral[:, 0] = 1e-5 * np.exp(-1j * WAVENUMBER * 8.0)
    dihedral[:, 3] = -dihedral[:, 0]
    dihedral[:, 1:3] = 1e-8
    np.savez(
        tmp_path / "dihedral.npz",
        s=dihedral,
        frequencies_hz=FREQUENCIES,
        sweep_values=np.array([1.41, 1.91, 2.41]),
        expressions=np.asarray(EXPRESSIONS),
        metadata=json.dumps({"sweep_variable": "copolarSpacingLambda"}),
    )
    main([str(tmp_path / "dihedral.npz"), "--out", str(tmp_path)])
    assert "All 3 values of copolarSpacingLambda give identical data" in capsys.readouterr().out
    written = {p.name for p in (tmp_path / "dihedral").iterdir()}
    assert {
        "signatures_vs_sweep.png",
        "signatures_vs_frequency.png",
        "signatures_vs_sweep.csv",
        "range_profiles.png",
    } <= written


def test_bistatic_angle_is_the_angle_at_the_target():
    target = [4.0, 0.0, 0.0]
    assert bistatic_angle_deg([0, 0, 0], [0, 0, 0], target) == pytest.approx(0.0)
    rx = [0.0, 4.0 * np.tan(np.radians(3.0)), 0.0]
    assert bistatic_angle_deg([0, 0, 0], rx, target) == pytest.approx(3.0)
    angles = bistatic_angle_deg(np.zeros((2, 3)), [[0, 0, 0], rx], [target, target])
    np.testing.assert_allclose(angles, [0.0, 3.0])
