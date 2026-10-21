#!/usr/bin/env python3
"""
Module 3: Offline polarimetric MIMO radar processing (pure NumPy, no AEDT).

Works on SBR+ sweeps exported by `sbr_results` (.npz): every channel S(<Rx port>,<Tx port>) over
frequency, for every value of a swept design variable. `MimoPolarimetryProcessor` runs the usual
MIMO radar chain on them and returns general products (range profiles, beamformed maps of total
power, polarimetric channel powers and Pauli features); `CanonicalTargetProcessor` adds a
single target's scattering matrix over frequency and its target-agnostic signatures (channel
powers, co-pol ratio and phase difference, co- minus cross-pol, reciprocity, Pauli, Cloude-Pottier).
Studies such as `bistaticity_analysis` combine these products; they do no signal processing of
their own.

The chain, one method each:

(a) Complex S-parameters: the export holds them as complex values (AEDT's real and imaginary
    parts, exact); `complex_from_db_deg` rebuilds them from magnitude (dB) and phase (deg) data.
(b) `frequency_window`: taper over the frequency samples.
(c) `range_profiles`: IFFT of the windowed sweep, with its phase referred to the centre frequency.
(d) `spatial_weights`: optional taper over each polarimetric channel's virtual aperture.
(e) Directions (azimuth, elevation), given by the caller; the virtual channel positions come from
    the export's builder port table (or `positions`).
(f) `beamform`: b_q = sum_p w_p x_q[p] conj(a_p(az, el)) / sum_p w_p.
(g) `range_azimuth_maps`, `azimuth_elevation_maps` and the map quantities `total_power`,
    `channel_powers` and `pauli_powers`.

A static scene swept in frequency has no slow time, hence no Doppler axis: maps are range-azimuth
and azimuth-elevation. `raw` forwards the exported data untouched.

Conventions
-----------
* Scattering matrices have shape (..., 2, 2) and are indexed S[..., rx, tx] with 0 = H and 1 = V,
  so S = [[S_HH, S_HV], [S_VH, S_VV]] and S_HV is the H-receive, V-transmit channel (AEDT's
  S(RxH,TxV)). The polarimetric channel q of S(Rx p, Tx q') is the string p + q', e.g. 'HV'.
* The array frame is the builder's global frame (= SBR+ AntennaCS): boresight +z, x horizontal
  (H), y vertical (V). A direction (az, el) is r = (cos el sin az, sin el, cos el cos az).
* AEDT's time convention is e^{+j w t}, so a path L contributes exp(-j k L). A Tx at p_t and an
  Rx at p_r see a far-field target in direction r with exp(+j k v . r) times common terms, where
  v = p_t + p_r is the virtual channel position. Hence the steering vector a_p = exp(+j k0 v_p . r)
  and the matched beamformer uses its conjugate.

Channels are read from their port names, both the hand-built style ('HornCluster1_RxH') and the
array builder's ('Linked1_Port_Rx_1_H'). When the export carries the builder's port table
(`sbr_results` attaches it), the port field signs are removed and the port phase centres give the
virtual channel positions.
"""

import json

import numpy as np
from scipy.signal import get_window

from pomaar.simulator.array_templates import (
    parse_port,
    polarimetric_quads,
    s_parameter_ports,
    select_quad,
)

SPEED_OF_LIGHT = 299792458.0

# Reference scattering matrices (BSA, S[rx, tx]) for fidelity metrics
REFERENCE_SPHERE = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=complex)
REFERENCE_TRIHEDRAL = REFERENCE_SPHERE  # odd bounce, like the sphere
REFERENCE_DIHEDRAL = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

POLARISATION_INDEX = {"H": 0, "V": 1}
POLARIMETRIC_CHANNELS = ("HH", "HV", "VH", "VV")


# ----------------------------------------------------------------------------------------------
# (a) Complex S-parameters and channel bookkeeping
# ----------------------------------------------------------------------------------------------


def complex_from_db_deg(magnitude_db, phase_deg):
    """S = 10^(M/20) exp(j pi/180 Phi) from a magnitude in dB and a phase in degrees."""
    return 10 ** (np.asarray(magnitude_db) / 20) * np.exp(1j * np.radians(phase_deg))


def assemble_scattering_matrix(s_hh, s_hv, s_vh, s_vv):
    """Stacks four channel arrays of equal shape into scattering matrices of shape (..., 2, 2)."""
    return np.stack([np.stack([s_hh, s_hv], axis=-1), np.stack([s_vh, s_vv], axis=-1)], axis=-2)


def channel_ports(expressions):
    """(receive port, transmit port) of every 'S(<receive>,<transmit>)' expression."""
    ports = []
    for expression in expressions:
        pair = s_parameter_ports(expression)
        if pair is None:
            raise ValueError(f"Expected S(receive, transmit): {expression}")
        ports.append(pair)
    return ports


def channel_layout(expressions):
    """Maps AEDT 'S(<Rx port>,<Tx port>)' expressions to (rx, tx) polarization index pairs."""
    layout = []
    for rx_port, tx_port in channel_ports(expressions):
        rx, tx = parse_port(rx_port), parse_port(tx_port)
        if not rx or not tx or rx[0] != "Rx" or tx[0] != "Tx":
            raise ValueError(
                f"Expected S(receive, transmit) with Rx/Tx H/V port names: S({rx_port},{tx_port})"
            )
        layout.append((POLARISATION_INDEX[rx[1]], POLARISATION_INDEX[tx[1]]))
    return layout


def quad_channel_indices(expressions, quad=None):
    """
    Indices of the S_HH, S_HV, S_VH, S_VV expressions of one polarimetric quad, as a 2x2 array
    [rx, tx], plus the quad's name and ports (see `array_templates.select_quad`; `quad` is None
    for the only quad, a quad name such as 'Rx_1/Tx_2', or a dict of the four ports).
    """
    pairs = channel_ports(expressions)
    ports = [port for pair in pairs for port in pair]
    name, ports_of = select_quad(ports, quad)
    index = {pair: i for i, pair in enumerate(pairs)}
    indices = np.zeros((2, 2), dtype=int)
    for rx_pol, rx in POLARISATION_INDEX.items():
        for tx_pol, tx in POLARISATION_INDEX.items():
            pair = (ports_of[f"Rx{rx_pol}"], ports_of[f"Tx{tx_pol}"])
            if pair not in index:
                raise ValueError(
                    f"Quad '{name}' lacks the {rx_pol}{tx_pol} channel S({pair[0]},{pair[1]})."
                )
            indices[rx, tx] = index[pair]
    return indices, name, ports_of


def to_scattering_matrices(s, expressions, quad=None):
    """
    (n_var, n_expr, n_freq) channel array -> (n_var, n_freq, 2, 2) matrices S[rx, tx] of one
    polarimetric quad (by default the only one in the export).
    """
    indices, _, _ = quad_channel_indices(expressions, quad)
    return np.moveaxis(np.asarray(s)[:, indices, :], -1, 1)


def port_sign_matrix(quad_ports, signs):
    """
    s_rx * s_tx per channel [rx, tx] from the per-port field signs (+1 along +x for H, +y for V).
    The measured matrix is diag(s_RxH, s_RxV) S diag(s_TxH, s_TxV), so multiplying it by this
    matrix (its own inverse) restores S.
    """
    rx = np.array([signs[quad_ports["RxH"]], signs[quad_ports["RxV"]]], dtype=float)
    tx = np.array([signs[quad_ports["TxH"]], signs[quad_ports["TxV"]]], dtype=float)
    return np.outer(rx, tx)


class IncompleteSweepError(ValueError):
    """Raised when some variations of an export lack data at some frequencies (NaN or inf)."""


def check_complete_sweep(s, frequencies_hz, sweep_values, sweep_variable=None):
    """
    Checks that every variation of `s` (n_variations, n_channels, n_frequencies) holds finite
    data at every frequency, and raises `IncompleteSweepError` otherwise, naming each incomplete
    variation with how many frequencies it has and their step.

    AEDT returns all variations on the union of their frequency grids and marks missing samples
    as NaN, so a variation that was solved on another grid (e.g. an old sweep), or not at all,
    shows up here.
    """
    s = np.asarray(s)
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    valid = np.all(np.isfinite(s), axis=1)  # (n_variations, n_frequencies)
    incomplete = np.flatnonzero(~valid.all(axis=-1))
    if incomplete.size == 0:
        return
    name = sweep_variable or "variation"
    values = np.arange(s.shape[0]) if sweep_values is None else np.asarray(sweep_values)

    def describe(index):
        solved = frequencies_hz[valid[index]]
        if solved.size == 0:
            return "no data at all"
        steps = np.unique(np.round(np.diff(solved) / 1e6, 6)) if solved.size > 1 else []
        step = f"step {steps[0]:g} MHz" if len(steps) == 1 else "irregular steps"
        return f"{solved.size} of {len(frequencies_hz)} frequencies ({step})"

    lines = [f"  {name} = {values[i]:g}: {describe(i)}" for i in incomplete]
    complete = [f"{values[i]:g}" for i in range(s.shape[0]) if i not in incomplete]
    raise IncompleteSweepError(
        f"{len(incomplete)} of {s.shape[0]} values of {name} lack data at some of the "
        f"{len(frequencies_hz)} frequencies, i.e. they were solved on another frequency grid or "
        "not at all:\n"
        + "\n".join(lines)
        + (f"\n  complete: {', '.join(complete)}" if complete else "")
        + "\nRe-solve every variation (the parametric sweep, not only the nominal setup) and "
        "export again."
    )


def identical_variations(s, sweep_values, sweep_variable=None, rtol=1e-9):
    """
    A warning naming the variations of `s` (n_variations, n_channels, n_frequencies) whose data
    equal those of another variation (to `rtol` of the export's largest magnitude), or None.

    Different values of a swept variable that give the same solution mean the variable changed
    nothing in the solved geometry, e.g. an SBR+ variable that is not mapped to the linked
    antenna. The data stay usable, so this only warns.
    """
    s = np.asarray(s)
    if s.shape[0] < 2:
        return None
    tolerance = rtol * np.max(np.abs(s))
    name = sweep_variable or "variation"
    values = np.arange(s.shape[0]) if sweep_values is None else np.asarray(sweep_values)
    groups, assigned = [], set()
    for i in range(s.shape[0]):
        if i in assigned:
            continue
        group = [i] + [
            j
            for j in range(i + 1, s.shape[0])
            if j not in assigned and np.max(np.abs(s[j] - s[i])) <= tolerance
        ]
        if len(group) > 1:
            groups.append(group)
            assigned.update(group)
    if not groups:
        return None
    if len(groups) == 1 and len(groups[0]) == s.shape[0]:
        return (
            f"All {s.shape[0]} values of {name} give identical data: the sweep changed nothing "
            f"in the solution (is {name} mapped to the linked antenna?)."
        )
    listed = "; ".join(", ".join(f"{values[i]:g}" for i in group) for group in groups)
    return (
        f"Some values of {name} give identical data ({listed}): those variations changed "
        "nothing in the solution."
    )


# ----------------------------------------------------------------------------------------------
# (b)-(f) Radar processing
# ----------------------------------------------------------------------------------------------


def frequency_window(n_frequencies, window="hann"):
    """Symmetric taper over the frequency samples; None or 'rect' is uniform."""
    if window in (None, "rect", "rectangular", "boxcar"):
        return np.ones(n_frequencies)
    return get_window(window, n_frequencies, fftbins=False)


def range_axis(frequencies_hz, n_range):
    """Range bins (m) of an `n_range`-point IFFT: c m / (2 n_range df), up to c / (2 df)."""
    df = frequencies_hz[1] - frequencies_hz[0]
    return SPEED_OF_LIGHT * np.arange(n_range) / (2 * n_range * df)


def range_compress(s, frequencies_hz, window="hann", oversampling=4):
    """
    IFFT of the windowed sweep along the last axis -> (range_m, profiles).

    A point scatterer of frequency-domain amplitude A at range R peaks at |A| in bin R. The phase
    of the profile refers to the centre frequency (the IFFT runs over the frequency index counted
    from the centre), which is what the steering vector's k0 assumes. Targets beyond the
    unambiguous range c / (2 df) fold back.
    """
    s = np.asarray(s)
    n_frequencies = s.shape[-1]
    n_range = int(oversampling * n_frequencies)
    taper = frequency_window(n_frequencies, window)
    profile = np.fft.ifft(s * taper, n=n_range, axis=-1) * n_range / taper.sum()
    centre_index = (n_frequencies - 1) / 2
    profile = profile * np.exp(-2j * np.pi * centre_index * np.arange(n_range) / n_range)
    return range_axis(frequencies_hz, n_range), profile


def direction_vectors(azimuth_deg, elevation_deg):
    """Unit vectors (n_az, n_el, 3) of the direction grid; boresight +z, azimuth towards +x."""
    az = np.radians(np.atleast_1d(np.asarray(azimuth_deg, dtype=float)))[:, None]
    el = np.radians(np.atleast_1d(np.asarray(elevation_deg, dtype=float)))[None, :]
    return np.stack(
        np.broadcast_arrays(np.cos(el) * np.sin(az), np.sin(el), np.cos(el) * np.cos(az)), axis=-1
    )


def steering_vectors(virtual_positions, azimuth_deg, elevation_deg, wavenumber):
    """
    a_p(az, el) = exp(+j k v_p . r) for virtual positions (..., n_p, 3) in metres, as an array
    (..., n_p, n_az, n_el) for a scalar k, or (..., n_p, n_k, n_az, n_el) for an array of k.
    """
    projection = np.einsum(
        "...pi,aei->...pae", virtual_positions, direction_vectors(azimuth_deg, elevation_deg)
    )
    k = np.asarray(wavenumber, dtype=float)
    if k.ndim == 0:
        return np.exp(1j * k * projection)
    return np.exp(1j * k[:, None, None] * projection[..., None, :, :])


def spatial_window(virtual_positions, window=None):
    """
    Weights (..., n_p) of a separable taper over the virtual aperture: the window function of
    each in-plane coordinate (x, y), evaluated over the aperture extent widened by one mean
    element spacing so that the edge channels keep a non-zero weight. None is uniform; 'hann',
    'hamming' and 'blackman' are available for arbitrary (sparse, 2D) virtual arrays.
    """
    positions = np.asarray(virtual_positions, dtype=float)
    weights = np.ones(positions.shape[:-1])
    if window is None:
        return weights
    shapes = {
        "hann": lambda u: 0.5 + 0.5 * np.cos(2 * np.pi * u),
        "hamming": lambda u: 0.54 + 0.46 * np.cos(2 * np.pi * u),
        "blackman": lambda u: 0.42 + 0.5 * np.cos(2 * np.pi * u) + 0.08 * np.cos(4 * np.pi * u),
    }
    if window not in shapes:
        raise ValueError(f"Unknown spatial window '{window}' (use {', '.join(shapes)} or None).")
    for axis in (0, 1):
        coordinate = positions[..., axis]
        low = coordinate.min(axis=-1, keepdims=True)
        high = coordinate.max(axis=-1, keepdims=True)
        extent = high - low
        n_unique = max(len(np.unique(np.round(coordinate, 9))), 2)
        width = np.where(extent > 0, extent * n_unique / (n_unique - 1), 1.0)
        u = (coordinate - (low + high) / 2) / width
        weights = weights * np.where(extent > 0, shapes[window](u), 1.0)
    return weights


def beamform(values, steering, weights):
    """
    b(az, el) = sum_p w_p x[p] conj(a_p(az, el)) / sum_p w_p, with values (n_var, n_p, n_x),
    steering (n_var, n_p, n_az, n_el) or (n_var, n_p, n_x, n_az, n_el) and weights (n_var, n_p).
    Returns (n_var, n_x, n_az, n_el).
    """
    weighted = values * weights[..., None] / weights.sum(axis=-1)[:, None, None]
    if steering.ndim == 4:
        return np.einsum("vpx,vpae->vxae", weighted, np.conj(steering))
    return np.einsum("vpx,vpxae->vxae", weighted, np.conj(steering))


# ----------------------------------------------------------------------------------------------
# (g) Map quantities
# ----------------------------------------------------------------------------------------------


def total_power(channel_maps):
    """Mean power over the polarimetric channels, mean_q |b_q|^2 (span / 4)."""
    return np.mean([np.abs(channel_maps[q]) ** 2 for q in channel_maps], axis=0)


def channel_powers(channel_maps, reference=None):
    """
    |b_q|^2 of every polarimetric channel, divided by one common reference so that their levels
    stay comparable: by default the peak of the total power (per variation, i.e. per leading
    index of the maps).
    """
    if reference is None:
        total = total_power(channel_maps)
        reference = total.reshape(total.shape[0], -1).max(axis=-1)
        reference = reference.reshape((-1,) + (1,) * (total.ndim - 1))
    return {q: np.abs(values) ** 2 / reference for q, values in channel_maps.items()}


def pauli_powers(channel_maps):
    """
    |k_i|^2 of the bistatic Pauli vector k = [HH+VV, HH-VV, HV+VH, j(HV-VH)] / sqrt(2) per pixel,
    as an array (..., 4); needs all four polarimetric channels.
    """
    matrices = assemble_scattering_matrix(*(channel_maps[q] for q in POLARIMETRIC_CHANNELS))
    return np.abs(pauli_vector(matrices)) ** 2


def db10(x):
    return 10 * np.log10(np.abs(x) + 1e-300)


def db20(x):
    return 20 * np.log10(np.abs(x) + 1e-300)


# ----------------------------------------------------------------------------------------------
# Per-matrix polarimetric signatures
# ----------------------------------------------------------------------------------------------


def bistatic_path(tx_position, rx_position, target_position):
    """Two-way spherical-wave path |p_tx - x| + |p_rx - x| (positions in metres, last axis xyz)."""
    tx_position, rx_position = np.asarray(tx_position), np.asarray(rx_position)
    target_position = np.asarray(target_position)
    return np.linalg.norm(tx_position - target_position, axis=-1) + np.linalg.norm(
        rx_position - target_position, axis=-1
    )


def bistatic_angle_deg(tx_position, rx_position, target_position):
    """Angle (deg) at the target between its directions to Tx and Rx; 0 when monostatic."""
    to_tx = np.asarray(tx_position, dtype=float) - np.asarray(target_position, dtype=float)
    to_rx = np.asarray(rx_position, dtype=float) - np.asarray(target_position, dtype=float)
    cosine = np.sum(to_tx * to_rx, axis=-1) / (
        np.linalg.norm(to_tx, axis=-1) * np.linalg.norm(to_rx, axis=-1)
    )
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def plane_wave_phase(frequencies_hz, tx_position, rx_position, direction):
    """
    Linear (plane-wave) MIMO phase term exp(+j k (p_tx + p_rx) . r_hat) of one Tx/Rx pair: the
    steering phase of its virtual channel for a far-field target in unit direction r_hat.
    """
    direction = np.asarray(direction, dtype=float)
    direction = direction / np.linalg.norm(direction)
    projection = (np.asarray(tx_position) + np.asarray(rx_position)) @ direction
    k = 2 * np.pi * np.asarray(frequencies_hz) / SPEED_OF_LIGHT
    return np.exp(1j * np.multiply.outer(projection, k))


def spherical_wave_phase(frequencies_hz, tx_position, rx_position, target_position):
    """
    Exact geometric phase exp(-j k (|p_tx - x| + |p_rx - x|)) of one Tx/Rx pair and focal point,
    the near-field counterpart of `plane_wave_phase`.
    """
    path = bistatic_path(tx_position, rx_position, target_position)
    k = 2 * np.pi * np.asarray(frequencies_hz) / SPEED_OF_LIGHT
    return np.exp(-1j * np.multiply.outer(path, k))


def copolar_ratio(scattering_matrix):
    """Co-polar ratio S_VV / S_HH: 0 dB, 0 deg for a sphere; 0 dB, 180 deg for a dihedral."""
    return scattering_matrix[..., 1, 1] / scattering_matrix[..., 0, 0]


def crosspolar_isolation(scattering_matrix):
    """
    Cross-polar leakage per transmit polarization, (S_VH / S_HH, S_HV / S_VV): the power that a
    transmitted H (V) wave puts into the V (H) receiver relative to the H (V) receiver.
    """
    s = scattering_matrix
    return s[..., 1, 0] / s[..., 0, 0], s[..., 0, 1] / s[..., 1, 1]


def reciprocity_ratio(scattering_matrix):
    """S_HV / S_VH, which a monostatic reciprocal measurement forces to 1 (0 dB, 0 deg)."""
    return scattering_matrix[..., 0, 1] / scattering_matrix[..., 1, 0]


def span(scattering_matrix):
    """Total power |S_HH|^2 + |S_HV|^2 + |S_VH|^2 + |S_VV|^2 (squared Frobenius norm)."""
    return np.sum(np.abs(scattering_matrix) ** 2, axis=(-2, -1))


def normalized_asymmetry(scattering_matrix):
    """|S_HV - S_VH| / ||S||_F, a reciprocity measure that stays finite when S_HV, S_VH -> 0."""
    s = scattering_matrix
    return np.abs(s[..., 0, 1] - s[..., 1, 0]) / np.sqrt(span(s))


def pauli_vector(scattering_matrix, bistatic=True):
    """
    Pauli target vector of shape (..., 4), or (..., 3) with bistatic=False.

    k = [S_HH + S_VV, S_HH - S_VV, S_HV + S_VH, j (S_HV - S_VH)] / sqrt(2), whose squared norm
    is the span. The fourth (antisymmetric) component vanishes for a monostatic reciprocal
    measurement. With bistatic=False it is dropped and the third becomes sqrt(2) times the mean
    cross-pol term, i.e. the monostatic k_p = [S_HH + S_VV, S_HH - S_VV, 2 S_HV] / sqrt(2) with
    S_HV symmetrized.
    """
    s = scattering_matrix
    s_hh, s_hv, s_vh, s_vv = s[..., 0, 0], s[..., 0, 1], s[..., 1, 0], s[..., 1, 1]
    components = [s_hh + s_vv, s_hh - s_vv, s_hv + s_vh]
    if bistatic:
        components.append(1j * (s_hv - s_vh))
    return np.stack(components, axis=-1) / np.sqrt(2.0)


def pauli_power_fractions(scattering_matrix):
    """Fractions of the span in the four bistatic Pauli components (they sum to 1)."""
    power = np.abs(pauli_vector(scattering_matrix, bistatic=True)) ** 2
    return power / power.sum(axis=-1, keepdims=True)


def polarimetric_fidelity(scattering_matrix, reference_matrix):
    """
    Phase- and scale-invariant similarity |<k, k_ref>| / (||k|| ||k_ref||) in [0, 1]; 1 when S
    is a complex multiple of the reference.
    """
    k = pauli_vector(scattering_matrix, bistatic=True)
    k_ref = pauli_vector(np.asarray(reference_matrix, dtype=complex), bistatic=True)
    inner = np.abs(np.sum(k * np.conj(k_ref), axis=-1))
    return inner / (np.linalg.norm(k, axis=-1) * np.linalg.norm(k_ref))


def alpha_angle(scattering_matrix):
    """Single-target Cloude alpha in degrees, cos(alpha) = |k_1| / ||k|| (0 sphere, 90 dihedral)."""
    k = pauli_vector(scattering_matrix, bistatic=True)
    ratio = np.abs(k[..., 0]) / np.linalg.norm(k, axis=-1)
    return np.degrees(np.arccos(np.clip(ratio, 0.0, 1.0)))


def coherency_matrix(target_vectors, axis=-2):
    """Averaged coherency matrix <k k^H> over `axis` of a stack of target vectors (..., N, d)."""
    k = np.moveaxis(target_vectors, axis, -2)
    return np.einsum("...ni,...nj->...ij", k, np.conj(k)) / k.shape[-2]


def cloude_pottier(coherency):
    """
    Cloude-Pottier parameters of (..., d, d) coherency matrices (bistatic 4x4 or monostatic 3x3):
    entropy H (log base d), anisotropy A = (l2 - l3) / (l2 + l3) of the eigenvalues sorted
    l1 >= l2 >= l3 >= ..., and the mean alpha (deg). Returns (H, A, mean_alpha). For one
    deterministic target averaged over frequency, H > 0 only if its polarimetric signature
    changes across the band.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(coherency)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    probabilities = eigenvalues / eigenvalues.sum(axis=-1, keepdims=True)
    dimension = coherency.shape[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(probabilities > 0, np.log(probabilities), 0.0)
        descending = np.sort(eigenvalues, axis=-1)[..., ::-1]
        anisotropy = np.nan_to_num(
            (descending[..., 1] - descending[..., 2]) / (descending[..., 1] + descending[..., 2])
        )
    entropy = -np.sum(probabilities * log_p, axis=-1) / np.log(dimension)
    alphas = np.arccos(np.clip(np.abs(eigenvectors[..., 0, :]), 0.0, 1.0))
    return entropy, anisotropy, np.degrees(np.sum(probabilities * alphas, axis=-1))


def reference_dihedral(roll_deg=0.0):
    """Dihedral rolled by `roll_deg` about the line of sight: R(roll) diag(1, -1) R(roll)^T."""
    c, s = np.cos(np.radians(roll_deg)), np.sin(np.radians(roll_deg))
    rotation = np.array([[c, -s], [s, c]])
    return (rotation @ REFERENCE_DIHEDRAL.real @ rotation.T).astype(complex)


def calibrate_with_sphere(matrices, sphere_matrices):
    """
    Complex channel-imbalance calibration of `matrices` from a sphere measured with the same
    channels (both (..., 2, 2)).

    A sphere's ideal matrix is the identity, so the co-polar gains are g_HH = S_HH and
    g_VV = S_VV of the sphere measurement. The cross-polar gains cannot be observed on a sphere;
    with g_pq = r_p t_q only their product g_HV g_VH = g_HH g_VV is known, so both are set to
    sqrt(g_HH g_VV). That is exact for a monostatic reciprocal system and an assumption for
    separate Tx/Rx antennas; a depolarizing reference (e.g. a 45 deg dihedral) is needed to
    resolve it. Cross-polar leakage is not corrected.
    """
    g_hh = sphere_matrices[..., 0, 0]
    g_vv = sphere_matrices[..., 1, 1]
    g_cross = np.sqrt(g_hh * g_vv)
    return matrices / assemble_scattering_matrix(g_hh, g_cross, g_cross, g_vv)


# ----------------------------------------------------------------------------------------------
# Processors
# ----------------------------------------------------------------------------------------------


class MimoPolarimetryProcessor:
    """
    MIMO radar processing of one exported SBR+ sweep: channels S(Rx, Tx) of shape
    (n_variations, n_channels, n_frequencies) over a uniform frequency grid, grouped into the
    polarimetric channels HH, HV, VH, VV, each with its virtual channels p (Tx/Rx pairs).

    Typical use:

        processor = MimoPolarimetryProcessor.from_npz("sphere_scattered.npz")
        range_m, profiles = processor.range_profiles()                  # (b), (c)
        maps = processor.azimuth_elevation_maps(np.arange(-60, 61), np.arange(-30, 31))
        power = total_power(maps["channels"])                          # (g)

    It expects the scattered (target-only) solution, the default of `sbr_results`.
    """

    def __init__(
        self,
        s,
        frequencies_hz,
        expressions,
        sweep_values=None,
        metadata=None,
        positions=None,
        correct_signs=True,
    ):
        """
        `s` is (n_variations, n_expressions, n_frequencies) complex and `expressions` its AEDT
        names S(<Rx port>,<Tx port>). `positions` maps port (excitation) names, or their
        suffixes such as 'TxH', to phase-centre positions in metres, (3,) or per variation
        (n_variations, 3); without it the export's builder port table is used, and channels
        without either get position 0 (no steering phase, exact only at boresight).
        """
        self.s = np.asarray(s, dtype=complex)
        self.frequencies_hz = np.asarray(frequencies_hz, dtype=float)
        self.expressions = [str(e) for e in expressions]
        n_variations = self.s.shape[0]
        self.sweep_values = (
            np.arange(n_variations, dtype=float)
            if sweep_values is None
            else np.asarray(sweep_values, dtype=float)
        )
        self.metadata = metadata or {}
        check_complete_sweep(
            self.s, self.frequencies_hz, self.sweep_values, self.metadata.get("sweep_variable")
        )
        self._raw = {
            "s": self.s.copy(),
            "frequencies_hz": self.frequencies_hz.copy(),
            "sweep_values": self.sweep_values.copy(),
            "expressions": list(self.expressions),
            "metadata": self.metadata,
        }

        self.channel_ports = channel_ports(self.expressions)
        self.channel_names = []
        for rx_port, tx_port in self.channel_ports:
            rx, tx = parse_port(rx_port), parse_port(tx_port)
            if not rx or not tx or rx[0] != "Rx" or tx[0] != "Tx":
                raise ValueError(f"Expected S(<Rx port>,<Tx port>): S({rx_port},{tx_port})")
            self.channel_names.append(rx[1] + tx[1])

        self.signs_corrected = False
        if correct_signs:
            self._correct_port_signs()
        self.positions_known = True
        self._port_positions = self._resolve_positions(positions)

    # --- loading and raw data -----------------------------------------------------------------

    @classmethod
    def from_npz(cls, path, quad=None, positions=None, correct_signs=True):
        """
        Loads an `sbr_results` export. `quad` (a name from `available_quads`, or a dict of the
        four ports) keeps only the four channels of one polarimetric quad; by default every
        channel is kept.
        """
        data = np.load(path)
        s = data["s"]
        expressions = list(data["expressions"])
        if quad is not None:
            indices, _, _ = quad_channel_indices(expressions, quad)
            keep = sorted(indices.ravel())
            s = s[:, keep]
            expressions = [expressions[i] for i in keep]
        return cls(
            s,
            data["frequencies_hz"],
            expressions,
            data["sweep_values"],
            json.loads(str(data["metadata"])),
            positions=positions,
            correct_signs=correct_signs,
        )

    @classmethod
    def available_quads(cls, path):
        """Names of the polarimetric quads in an `sbr_results` export."""
        expressions = np.load(path)["expressions"]
        return list(polarimetric_quads([p for pair in channel_ports(expressions) for p in pair]))

    def raw(self):
        """The exported data as loaded, untouched."""
        return dict(self._raw)

    @property
    def port_table(self):
        """Builder port table rows of the export ({excitation: row}); empty without one."""
        rows = self.metadata.get("ports") or []
        return {row["excitation"]: row for row in rows if row.get("excitation")}

    def _correct_port_signs(self):
        """
        Takes the builder's port field signs out of every channel (S_measured = s_rx s_tx S),
        when the export's port table covers all ports.
        """
        table = self.port_table
        ports = {port for pair in self.channel_ports for port in pair}
        if not ports or not all(table.get(port, {}).get("sign") for port in ports):
            return
        signs = np.array([table[rx]["sign"] * table[tx]["sign"] for rx, tx in self.channel_ports])
        self.s = self.s * signs[None, :, None]
        self.signs_corrected = True

    def _resolve_positions(self, positions):
        """{port: (n_variations, 3) metres} for every port of the export."""
        n_variations = self.s.shape[0]
        table = self.port_table
        resolved = {}
        for port in {port for pair in self.channel_ports for port in pair}:
            value = None
            if positions:
                key = (
                    port
                    if port in positions
                    else next(
                        (k for k in positions if port.endswith(f"_{k}") or port.endswith(k)), None
                    )
                )
                value = positions.get(key) if key is not None else None
            if value is None and port in table and table[port].get("phase_centre_mm"):
                value = np.asarray(table[port]["phase_centre_mm"], dtype=float) / 1000.0
            if value is None:
                self.positions_known = False
                value = np.zeros(3)
            resolved[port] = np.broadcast_to(np.asarray(value, dtype=float), (n_variations, 3))
        return resolved

    # --- channels and geometry ----------------------------------------------------------------

    @property
    def centre_frequency_hz(self):
        return 0.5 * (self.frequencies_hz[0] + self.frequencies_hz[-1])

    @property
    def centre_wavenumber(self):
        return 2 * np.pi * self.centre_frequency_hz / SPEED_OF_LIGHT

    @property
    def polarimetric_channels(self):
        """{q: [channel indices]} for the polarimetric channels present (HH, HV, VH, VV)."""
        groups = {}
        for index, name in enumerate(self.channel_names):
            groups.setdefault(name, []).append(index)
        return {q: groups[q] for q in POLARIMETRIC_CHANNELS if q in groups}

    def virtual_positions(self, q):
        """Virtual channel positions v_p = p_tx + p_rx (n_variations, n_p, 3) of channel q."""
        return np.stack(
            [
                self._port_positions[self.channel_ports[i][1]]
                + self._port_positions[self.channel_ports[i][0]]
                for i in self.polarimetric_channels[q]
            ],
            axis=1,
        )

    def spatial_weights(self, q, window=None):
        """(d) Taper weights (n_variations, n_p) over the virtual aperture of channel q."""
        return spatial_window(self.virtual_positions(q), window)

    # --- (b), (c) range compression -----------------------------------------------------------

    def frequency_window(self, window="hann"):
        """(b) The taper applied over the frequency samples."""
        return frequency_window(len(self.frequencies_hz), window)

    def range_profiles(self, window="hann", oversampling=4):
        """
        (c) Range profiles of every channel: (range_m, profiles (n_variations, n_channels,
        n_range)), with the phase referred to the centre frequency.
        """
        return range_compress(self.s, self.frequencies_hz, window, oversampling)

    # --- (e), (f) beamforming -----------------------------------------------------------------

    def beamform(self, values, azimuth_deg, elevation_deg, spatial_window=None, wavenumber=None):
        """
        (f) Beamforms channel data `values` (n_variations, n_channels, n_x) (range profiles, or
        the frequency samples with `wavenumber` = their k) towards every (azimuth, elevation),
        separately per polarimetric channel: {q: (n_variations, n_x, n_az, n_el)}.
        """
        k = self.centre_wavenumber if wavenumber is None else wavenumber
        result = {}
        for q, indices in self.polarimetric_channels.items():
            steering = steering_vectors(self.virtual_positions(q), azimuth_deg, elevation_deg, k)
            weights = self.spatial_weights(q, spatial_window)
            result[q] = beamform(values[:, indices], steering, weights)
        return result

    # --- (g) maps ----------------------------------------------------------------------------

    def range_azimuth_maps(
        self,
        azimuth_deg,
        elevation_deg=0.0,
        window="hann",
        spatial_window=None,
        oversampling=4,
    ):
        """
        Range-azimuth maps at one elevation: {'range_m', 'azimuth_deg', 'channels': {q:
        complex (n_variations, n_range, n_az)}}. Feed 'channels' to `total_power`,
        `channel_powers` or `pauli_powers`.
        """
        range_m, profiles = self.range_profiles(window, oversampling)
        beams = self.beamform(profiles, azimuth_deg, [elevation_deg], spatial_window)
        return {
            "range_m": range_m,
            "azimuth_deg": np.atleast_1d(azimuth_deg),
            "elevation_deg": float(elevation_deg),
            "channels": {q: b[..., 0] for q, b in beams.items()},
        }

    def peak_range_bins(self, profiles):
        """Range bin (per variation) of the strongest return, summed incoherently over channels."""
        return np.argmax(np.sum(np.abs(profiles) ** 2, axis=1), axis=-1)

    def azimuth_elevation_maps(
        self,
        azimuth_deg,
        elevation_deg,
        range_m=None,
        window="hann",
        spatial_window=None,
        oversampling=4,
    ):
        """
        Azimuth-elevation maps at one range: {'range_m' (per variation), 'azimuth_deg',
        'elevation_deg', 'channels': {q: complex (n_variations, n_az, n_el)}}. The range is
        `range_m`, or the strongest return of each variation (`peak_range_bins`).
        """
        ranges, profiles = self.range_profiles(window, oversampling)
        if range_m is None:
            bins = self.peak_range_bins(profiles)
        else:
            bins = np.full(profiles.shape[0], np.argmin(np.abs(ranges - range_m)))
        selected = profiles[np.arange(profiles.shape[0]), :, bins][..., None]
        beams = self.beamform(selected, azimuth_deg, elevation_deg, spatial_window)
        return {
            "range_m": ranges[bins],
            "azimuth_deg": np.atleast_1d(azimuth_deg),
            "elevation_deg": np.atleast_1d(elevation_deg),
            "channels": {q: b[:, 0] for q, b in beams.items()},
        }

    def inject_noise(self, signal_to_noise_ratio_db, rng=None):
        """Adds complex white Gaussian noise at the given SNR relative to the mean channel power."""
        rng = np.random.default_rng(rng)
        signal_power = np.mean(np.abs(self.s) ** 2)
        noise_std = np.sqrt(signal_power / 10 ** (signal_to_noise_ratio_db / 10) / 2)
        self.s = self.s + noise_std * (
            rng.standard_normal(self.s.shape) + 1j * rng.standard_normal(self.s.shape)
        )


class CanonicalTargetProcessor(MimoPolarimetryProcessor):
    """
    MIMO processing plus the polarimetric description of a single target (e.g. a canonical
    scatterer): its 2x2 scattering matrix over frequency, beamformed towards its direction, and
    target-agnostic signatures of that matrix. Nothing here assumes what the target is; a
    reference target (`fidelity`) or a classification (`bounce_parity`) is an explicit opt-in.

        processor = CanonicalTargetProcessor.from_npz("dihedral.npz")
        matrices = processor.target_matrices()                   # (n_var, n_freq, 2, 2)
        features = processor.signatures(matrices=matrices)
    """

    def target_matrices(self, azimuth_deg=0.0, elevation_deg=0.0, spatial_window=None):
        """
        The target's scattering matrices S[rx, tx] (n_variations, n_frequencies, 2, 2): every
        polarimetric channel beamformed towards (azimuth, elevation) in the frequency domain,
        with the exact k of every frequency. With one virtual channel per polarimetric channel
        (a cluster) this is the plane-wave compensation of each channel towards the target,
        and at boresight the channels as exported.
        """
        k = 2 * np.pi * self.frequencies_hz / SPEED_OF_LIGHT
        beams = self.beamform(self.s, [azimuth_deg], [elevation_deg], spatial_window, k)
        missing = [q for q in POLARIMETRIC_CHANNELS if q not in beams]
        if missing:
            raise ValueError(f"The export lacks the polarimetric channel(s) {', '.join(missing)}.")
        return assemble_scattering_matrix(*(beams[q][..., 0, 0] for q in POLARIMETRIC_CHANNELS))

    @staticmethod
    def copol_crosspol_difference_db(matrices):
        """
        Co-polar minus cross-polar return per transmit polarization in dB,
        (|S_HH| - |S_VH| for Tx H, |S_VV| - |S_HV| for Tx V).
        """
        cross_h, cross_v = crosspolar_isolation(matrices)
        return -db20(cross_h), -db20(cross_v)

    @staticmethod
    def copol_phase_difference_deg(matrices):
        """Co-polar phase difference arg(S_VV S_HH*) in degrees."""
        return np.angle(copolar_ratio(matrices), deg=True)

    @classmethod
    def bounce_parity(cls, matrices):
        """
        Opt-in classification: 'odd' where |arg(S_VV S_HH*)| < 90 deg (sphere, trihedral, plate),
        else 'even' (dihedral). It assumes a dominant single-, double- or triple-bounce target.
        """
        return np.where(np.abs(cls.copol_phase_difference_deg(matrices)) < 90.0, "odd", "even")

    @staticmethod
    def fidelity(matrices, reference_matrix):
        """Opt-in comparison with a chosen reference matrix (`polarimetric_fidelity`)."""
        return polarimetric_fidelity(matrices, reference_matrix)

    def signatures(self, azimuth_deg=0.0, elevation_deg=0.0, matrices=None):
        """
        Target-agnostic signatures of `matrices` (or of the target matrices towards
        (azimuth, elevation)), per variation and frequency unless noted:

        * `power_<pq>_db`: |S_pq|^2 of every channel, and `span_db`: the total |S|_F^2;
        * `copol_ratio_db`, `copol_phase_difference_deg`: |S_VV / S_HH|, arg(S_VV S_HH*);
        * `co_minus_cross_tx<h|v>_db`, `co_cross_phase_tx<h|v>_deg`: co- to cross-polar ratio
          per transmit polarization (the two receivers of one transmission, as in a weather
          radar's linear depolarization ratio), |S_co / S_cross| and arg(S_co S_cross*);
        * `reciprocity_db`, `reciprocity_phase_deg`, `asymmetry_db`: S_HV / S_VH and
          |S_HV - S_VH| / |S|_F;
        * `pauli_fractions_db`: the four bistatic Pauli power fractions;
        * `alpha_deg`: the Pauli angle of the single matrix, cos a = |k_1| / |k|;
        * `entropy`, `anisotropy`, `mean_alpha_deg` (per variation): Cloude-Pottier parameters
          of the bistatic 4x4 coherency matrix averaged over the band.
        """
        s = self.target_matrices(azimuth_deg, elevation_deg) if matrices is None else matrices
        co_minus_cross_h, co_minus_cross_v = self.copol_crosspol_difference_db(s)
        # arg(S_co S_cross*) = -arg(S_cross / S_co)
        cross_h, cross_v = crosspolar_isolation(s)
        entropy, anisotropy, mean_alpha = cloude_pottier(
            coherency_matrix(pauli_vector(s, bistatic=True))
        )
        ratio = copolar_ratio(s)
        reciprocity = reciprocity_ratio(s)
        powers = {
            f"power_{q.lower()}_db": db20(s[..., "HV".index(q[0]), "HV".index(q[1])])
            for q in POLARIMETRIC_CHANNELS
        }
        return {
            **powers,
            "span_db": db10(span(s)),
            "copol_ratio_db": db20(ratio),
            "copol_phase_difference_deg": np.angle(ratio, deg=True),
            "co_minus_cross_txh_db": co_minus_cross_h,
            "co_cross_phase_txh_deg": -np.angle(cross_h, deg=True),
            "co_minus_cross_txv_db": co_minus_cross_v,
            "co_cross_phase_txv_deg": -np.angle(cross_v, deg=True),
            "reciprocity_db": db20(reciprocity),
            "reciprocity_phase_deg": np.angle(reciprocity, deg=True),
            "asymmetry_db": db20(normalized_asymmetry(s)),
            "pauli_fractions_db": db10(pauli_power_fractions(s)),
            "alpha_deg": alpha_angle(s),
            "entropy": entropy,
            "anisotropy": anisotropy,
            "mean_alpha_deg": mean_alpha,
        }
