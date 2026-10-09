#!/usr/bin/env python3
"""
Module 3: Offline polarimetric processor (pure NumPy, no AEDT).

Works from SBR+ sweeps exported by `sbr_results` (.npz). The module-level functions are the
verified routines of the bistaticity study (see the README): range gating, geometric phase
compensation, and per-matrix polarimetric signatures and decompositions. `MimoPolarimetryProcessor`
chains them for one exported sweep.

Scattering matrices have shape (..., 2, 2) and are indexed S[..., rx, tx] with 0 = H and 1 = V,
so S = [[S_HH, S_HV], [S_VH, S_VV]] and S_HV is the H-receive, V-transmit channel (AEDT's
S(RxH,TxV)).
"""

import json
import re

import numpy as np
from scipy.signal import get_window

SPEED_OF_LIGHT = 299792458.0

# Reference scattering matrices (BSA, S[rx, tx]) for fidelity metrics
REFERENCE_SPHERE = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=complex)
REFERENCE_DIHEDRAL = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


POLARISATION_INDEX = {"H": 0, "V": 1}


def assemble_scattering_matrix(s_hh, s_hv, s_vh, s_vv):
    """Stacks four channel arrays of equal shape into scattering matrices of shape (..., 2, 2)."""
    return np.stack([np.stack([s_hh, s_hv], axis=-1), np.stack([s_vh, s_vv], axis=-1)], axis=-2)


def channel_layout(expressions):
    """Maps AEDT 'S(<..>Rx<p>..,<..>Tx<q>..)' expressions to (rx, tx) polarisation index pairs."""
    layout = []
    for expression in expressions:
        match = re.fullmatch(r"S\(\w*?Rx([HV])\w*,\w*?Tx([HV])\w*\)", str(expression))
        if not match:
            raise ValueError(f"Expected S(receive, transmit) with Rx/Tx port names: {expression}")
        layout.append((POLARISATION_INDEX[match[1]], POLARISATION_INDEX[match[2]]))
    return layout


def to_scattering_matrices(s, expressions):
    """(n_var, n_expr, n_freq) channel array -> (n_var, n_freq, 2, 2) matrices S[rx, tx]."""
    layout = channel_layout(expressions)
    if sorted(layout) != [(0, 0), (0, 1), (1, 0), (1, 1)]:
        raise ValueError(f"Expected the four HH/HV/VH/VV channels, got {list(expressions)}")
    matrices = np.zeros((s.shape[0], s.shape[2], 2, 2), dtype=complex)
    for expression_index, (rx, tx) in enumerate(layout):
        matrices[:, :, rx, tx] = s[:, expression_index, :]
    return matrices


def delay_profile(s_freq, frequencies_hz, window="blackmanharris", n_fft=None):
    """
    Transforms a uniformly swept response to the two-way path domain.

    Returns (path_m, profile) with path_m centred on zero (fftshift order) and profile scaled so
    that a single scatterer of frequency-domain amplitude A peaks at |A|. Paths are ambiguous
    modulo c/df, so a target farther than that wraps around.
    """
    n_freq = s_freq.shape[-1]
    n_fft = n_fft or n_freq
    taper = get_window(window, n_freq, fftbins=False)
    df = frequencies_hz[1] - frequencies_hz[0]
    profile = np.fft.ifft(s_freq * taper, n=n_fft, axis=-1) * n_fft / taper.sum()
    path_m = np.fft.fftfreq(n_fft, d=df) * SPEED_OF_LIGHT
    return np.fft.fftshift(path_m), np.fft.fftshift(profile, axes=-1)


def wrap_path(path_m, frequencies_hz):
    """Wraps an absolute two-way path into the centred unambiguous window of a frequency sweep."""
    unambiguous_m = SPEED_OF_LIGHT / (frequencies_hz[1] - frequencies_hz[0])
    return (path_m + unambiguous_m / 2) % unambiguous_m - unambiguous_m / 2


def time_gate(s_freq, frequencies_hz, centre_path_m, half_width_m, window="blackmanharris"):
    """
    Range-gates a swept response around one two-way path (VNA-style, gate-normalised).

    The sweep is tapered with `window`, transformed to the path domain, multiplied by a
    rectangular gate of +-`half_width_m` centred on `centre_path_m` (wrapped into the unambiguous
    window), and transformed back. The result is then divided by the same processing applied to
    an ideal point reflector exp(-j k centre_path_m) and multiplied by that reflector again. This
    gate normalisation removes the taper and the truncation of the target's own main lobe, so a
    point scatterer at `centre_path_m` comes out unbiased; only leakage of other paths into the
    gate remains. The band edges stay less reliable (see `gate_valid_band`).
    """
    n_freq = s_freq.shape[-1]
    taper = get_window(window, n_freq, fftbins=False)
    df = frequencies_hz[1] - frequencies_hz[0]
    path_m = np.fft.fftfreq(n_freq, d=df) * SPEED_OF_LIGHT
    unambiguous_m = SPEED_OF_LIGHT / df
    offset = (path_m - wrap_path(centre_path_m, frequencies_hz) + unambiguous_m / 2) % unambiguous_m
    gate = np.abs(offset - unambiguous_m / 2) <= half_width_m

    def apply_gate(response):
        return np.fft.fft(np.fft.ifft(response * taper, axis=-1) * gate, axis=-1)

    reflector = np.exp(-2j * np.pi * np.asarray(frequencies_hz) * centre_path_m / SPEED_OF_LIGHT)
    return apply_gate(s_freq) / apply_gate(reflector) * reflector


def gate_valid_band(frequencies_hz, window="blackmanharris", min_taper=0.5):
    """Boolean mask of the frequencies where the gating taper exceeds `min_taper` of its peak."""
    taper = get_window(window, len(frequencies_hz), fftbins=False)
    return taper >= min_taper * taper.max()


def delay_window_peak(
    s_freq, frequencies_hz, centre_path_m, offsets_m, window="blackmanharris", oversample=16
):
    """
    Peak |delay profile| within centre_path_m + [offsets_m[0], offsets_m[1]] (wrapped).

    With offsets (-half_width, +half_width) it is the gated target level; with an interval beside
    the gate (a guard band) it estimates what leaks into the gate from elsewhere, e.g. the skirt
    of a much stronger, dispersive direct-coupling term. Their ratio is the gate's
    signal-to-contamination ratio.
    """
    path_m, profile = delay_profile(
        s_freq, frequencies_hz, window=window, n_fft=oversample * s_freq.shape[-1]
    )
    offset = wrap_path(path_m - centre_path_m, frequencies_hz)
    inside = (offset >= offsets_m[0]) & (offset <= offsets_m[1])
    return np.max(np.abs(profile[..., inside]), axis=-1)


def bistatic_path(tx_position, rx_position, target_position):
    """Two-way spherical-wave path |p_tx - x| + |p_rx - x| (positions in metres, last axis xyz)."""
    tx_position, rx_position = np.asarray(tx_position), np.asarray(rx_position)
    target_position = np.asarray(target_position)
    return np.linalg.norm(tx_position - target_position, axis=-1) + np.linalg.norm(
        rx_position - target_position, axis=-1
    )


def plane_wave_phase(frequencies_hz, tx_position, rx_position, direction):
    """
    Linear (plane-wave) MIMO phase term exp(+j k (p_tx + p_rx) . r_hat) of one Tx/Rx pair.

    This is the far-field steering phase of the virtual element at (p_tx + p_rx) / 2, relative to
    the array origin, for a target in unit direction r_hat. Dividing a channel by it moves that
    channel's phase reference to the origin.
    """
    direction = np.asarray(direction, dtype=float)
    direction = direction / np.linalg.norm(direction)
    projection = (np.asarray(tx_position) + np.asarray(rx_position)) @ direction
    k = 2 * np.pi * np.asarray(frequencies_hz) / SPEED_OF_LIGHT
    return np.exp(1j * np.multiply.outer(projection, k))


def spherical_wave_phase(frequencies_hz, tx_position, rx_position, target_position):
    """
    Exact geometric phase exp(-j k (|p_tx - x| + |p_rx - x|)) of one Tx/Rx pair and focal point.

    Dividing a channel by it removes the linear term and the Fresnel (quadratic) term together,
    i.e. it focuses the channel on x. What remains is the scattering itself plus the antenna
    response in the direction of x.
    """
    path = bistatic_path(tx_position, rx_position, target_position)
    k = 2 * np.pi * np.asarray(frequencies_hz) / SPEED_OF_LIGHT
    return np.exp(-1j * np.multiply.outer(path, k))


def copolar_ratio(scattering_matrix):
    """Co-polar ratio S_VV / S_HH: 0 dB, 0 deg for a sphere; 0 dB, 180 deg for a dihedral."""
    return scattering_matrix[..., 1, 1] / scattering_matrix[..., 0, 0]


def crosspolar_isolation(scattering_matrix):
    """
    Cross-polar leakage per transmit polarisation, (S_VH / S_HH, S_HV / S_VV).

    The first ratio is the power that a transmitted H wave puts into the V receiver relative to
    the H receiver, the second the same for a transmitted V wave.
    """
    s = scattering_matrix
    return s[..., 1, 0] / s[..., 0, 0], s[..., 0, 1] / s[..., 1, 1]


def reciprocity_ratio(scattering_matrix):
    """S_HV / S_VH, which a monostatic reciprocal measurement forces to 1 (0 dB, 0 deg)."""
    return scattering_matrix[..., 0, 1] / scattering_matrix[..., 1, 0]


def span(scattering_matrix):
    """Total power |S_HH|^2 + |S_HV|^2 + |S_VH|^2 + |S_VV|^2 (squared Frobenius norm)."""
    return np.sum(np.abs(scattering_matrix) ** 2, axis=(-2, -1))


def normalised_asymmetry(scattering_matrix):
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
    S_HV symmetrised.
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
    Phase- and scale-invariant similarity |<k, k_ref>| / (||k|| ||k_ref||) in [0, 1].

    It is 1 when S is a complex multiple of the reference. 1 - fidelity is the infidelity; in dB
    it is the floor that a bistatic configuration imposes on telling the reference mechanism
    apart from others.
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


def entropy_alpha(coherency):
    """
    Cloude-Pottier entropy H (log base d) and mean alpha (deg) of (..., d, d) coherency matrices.

    Works for the monostatic 3x3 and the bistatic 4x4 matrices. Averaging the Pauli vectors of one
    deterministic target over frequency gives H > 0 only if its polarimetric signature changes
    across the band.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(coherency)
    eigenvalues = np.clip(eigenvalues, 0.0, None)
    probabilities = eigenvalues / eigenvalues.sum(axis=-1, keepdims=True)
    dimension = coherency.shape[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(probabilities > 0, np.log(probabilities), 0.0)
    entropy = -np.sum(probabilities * log_p, axis=-1) / np.log(dimension)
    alphas = np.arccos(np.clip(np.abs(eigenvectors[..., 0, :]), 0.0, 1.0))
    return entropy, np.degrees(np.sum(probabilities * alphas, axis=-1))


class MimoPolarimetryProcessor:
    """
    Polarimetric processing of one exported SBR+ sweep: scattering matrices of shape
    (n_variations, n_frequencies, 2, 2) over a uniform frequency grid.

    Typical use:

        processor = MimoPolarimetryProcessor.from_npz("sphere_scattered.npz")
        processor.compensate_geometry(positions_per_variation, focal_point)  # optional
        signatures = processor.signatures(reference=REFERENCE_SPHERE)

    Use the `Setup : Sweep_Scattered` export (target only). For total-field exports, pass the
    `Sweep_Incident` or no-target export as `background`, or range-gate with `gate`.
    """

    def __init__(self, scattering_matrices, frequencies_hz, sweep_values=None, metadata=None):
        self.scattering_matrices = np.asarray(scattering_matrices, dtype=complex)
        self.frequencies_hz = np.asarray(frequencies_hz, dtype=float)
        n_variations = self.scattering_matrices.shape[0]
        self.sweep_values = (
            np.arange(n_variations, dtype=float)
            if sweep_values is None
            else np.asarray(sweep_values)
        )
        self.metadata = metadata or {}
        self.valid_band = np.ones(len(self.frequencies_hz), dtype=bool)

    @classmethod
    def from_npz(cls, path, background_path=None):
        """Loads an `sbr_results` export, optionally subtracting a background export coherently."""
        data = np.load(path)
        s = data["s"]
        if background_path:
            background = np.load(background_path)
            for key in ("frequencies_hz", "sweep_values", "expressions"):
                if not np.array_equal(background[key], data[key]):
                    raise ValueError(f"Background sweep differs from the target sweep in '{key}'")
            s = s - background["s"]
        return cls(
            to_scattering_matrices(s, data["expressions"]),
            data["frequencies_hz"],
            data["sweep_values"],
            json.loads(str(data["metadata"])),
        )

    def gate(self, centre_paths_m, half_width_m, window=("kaiser", 14), min_taper=0.5):
        """
        Range-gates every channel around one two-way path per variation (`time_gate`) and limits
        `valid_band` to where the taper is at least `min_taper`. Only needed when the data still
        hold other paths, e.g. direct coupling in a total-field export or several scatterers.
        """
        centre_paths_m = np.broadcast_to(centre_paths_m, self.sweep_values.shape)
        for index, centre in enumerate(centre_paths_m):
            channels = np.moveaxis(self.scattering_matrices[index], 0, -1)
            gated = time_gate(channels, self.frequencies_hz, centre, half_width_m, window=window)
            self.scattering_matrices[index] = np.moveaxis(gated, -1, 0)
        self.valid_band &= gate_valid_band(self.frequencies_hz, window=window, min_taper=min_taper)

    def compensate_geometry(self, positions, focal_point):
        """
        Divides every channel by its spherical-wave phase to `focal_point` (`spherical_wave_phase`),
        removing the linear and Fresnel terms of the bistatic geometry.

        `positions` is one dict per variation with keys 'TxH', 'TxV', 'RxH', 'RxV' (metres), or a
        single dict used for all variations.
        """
        if isinstance(positions, dict):
            positions = [positions] * len(self.sweep_values)
        names = "HV"
        for index, element_positions in enumerate(positions):
            for rx in (0, 1):
                for tx in (0, 1):
                    phase = spherical_wave_phase(
                        self.frequencies_hz,
                        element_positions[f"Tx{names[tx]}"],
                        element_positions[f"Rx{names[rx]}"],
                        focal_point,
                    )
                    self.scattering_matrices[index, :, rx, tx] /= phase

    def calibrate_with_sphere(self, sphere):
        """
        Complex channel-imbalance calibration from a sphere measured with the same channels.

        A sphere's ideal matrix is the identity, so the co-polar gains are g_HH = S_HH and
        g_VV = S_VV of the sphere measurement (per variation and frequency). The cross-polar gains
        cannot be observed on a sphere; with g_pq = r_p t_q only their product g_HV g_VH = g_HH g_VV
        is known, so both are set to sqrt(g_HH g_VV). That is exact for a monostatic reciprocal
        system and an assumption for separate Tx/Rx antennas; a depolarising reference (e.g. a 45
        deg dihedral) is needed to resolve it. Cross-polar leakage is not corrected.
        """
        g_hh = sphere.scattering_matrices[..., 0, 0]
        g_vv = sphere.scattering_matrices[..., 1, 1]
        g_cross = np.sqrt(g_hh * g_vv)
        gains = assemble_scattering_matrix(g_hh, g_cross, g_cross, g_vv)
        self.scattering_matrices = self.scattering_matrices / gains

    def inject_noise(self, signal_to_noise_ratio_db, rng=None):
        """Adds complex white Gaussian noise at the given SNR relative to the mean channel power."""
        rng = np.random.default_rng(rng)
        signal_power = np.mean(np.abs(self.scattering_matrices) ** 2)
        noise_std = np.sqrt(signal_power / 10 ** (signal_to_noise_ratio_db / 10) / 2)
        shape = self.scattering_matrices.shape
        self.scattering_matrices = self.scattering_matrices + noise_std * (
            rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
        )

    def signatures(self, reference=REFERENCE_SPHERE):
        """
        All signatures of the README table, per variation and frequency (arrays of shape
        (n_variations, n_frequencies[, 4])), plus band-averaged entropy / mean alpha per variation
        over `valid_band`. Ratios are linear complex values; *_db keys are in dB.
        """
        s = self.scattering_matrices
        tx_h_leak, tx_v_leak = crosspolar_isolation(s)
        in_band = s[:, self.valid_band]
        entropy_bistatic, mean_alpha_bistatic = entropy_alpha(
            coherency_matrix(pauli_vector(in_band, bistatic=True))
        )
        entropy_monostatic, mean_alpha_monostatic = entropy_alpha(
            coherency_matrix(pauli_vector(in_band, bistatic=False))
        )
        return {
            "copol_ratio": copolar_ratio(s),
            "xpol_txh": tx_h_leak,
            "xpol_txv": tx_v_leak,
            "reciprocity": reciprocity_ratio(s),
            "asymmetry": normalised_asymmetry(s),
            "pauli_fractions": pauli_power_fractions(s),
            "infidelity": 1 - polarimetric_fidelity(s, reference),
            "alpha_deg": alpha_angle(s),
            "entropy_bistatic": entropy_bistatic,
            "mean_alpha_bistatic_deg": mean_alpha_bistatic,
            "entropy_monostatic": entropy_monostatic,
            "mean_alpha_monostatic_deg": mean_alpha_monostatic,
        }
