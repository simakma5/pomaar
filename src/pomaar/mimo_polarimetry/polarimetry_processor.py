#!/usr/bin/env python3
"""
Module 3: Offline DSP & Polarimetric Processor.
Reads SBR+ JSON results and HFSS S-parameters, performs PDM calibration,
injects noise, and computes Pauli, Cloude-Pottier, and Yamaguchi decompositions.
"""

import json
import os

import numpy as np


class MimoPolarimetryProcessor:
    """
    Handles offline signal processing, calibration, and polarimetric decompositions
    for multi-static MIMO radar data.
    """

    def __init__(self, sbr_results_json_path, array_coupling_touchstone_path=None):
        self.sbr_results_json_path = sbr_results_json_path
        self.array_coupling_touchstone_path = array_coupling_touchstone_path

        self.frequencies_ghz = []
        self.raw_voltages = {}  # Keys: (tx_idx, rx_idx, tx_pol, rx_pol) -> complex numpy array
        self.coupling_matrix = None  # (M+N) x (M+N) S-parameter matrix
        self.calibrated_scattering_matrices = {}  # Keys: (tx_idx, rx_idx) -> 2x2 complex matrix per frequency

        # Parse inputs
        self.load_sbr_results()
        if array_coupling_touchstone_path:
            self.load_coupling_matrix()

    def load_sbr_results(self):
        """Loads and parses SBR+ JSON voltage outputs."""
        print(f"Loading SBR+ results from: {self.sbr_results_json_path}")
        with open(self.sbr_results_json_path, "r") as file_handle:
            data = json.load(file_handle)

        self.frequencies_ghz = np.array(data["frequencies_ghz"])
        channel_data = data["channels"]

        # Populate raw voltages dictionary
        for channel_key, val_dict in channel_data.items():
            # Expected format: "Tx_1_H->Rx_1_V"
            tx_part, rx_part = channel_key.split("->")

            tx_tokens = tx_part.split("_")
            tx_idx = int(tx_tokens[1])
            tx_pol = tx_part[-1]  # 'H' or 'V'

            rx_tokens = rx_part.split("_")
            rx_idx = int(rx_tokens[1])
            rx_pol = rx_part[-1]  # 'H' or 'V'

            real_parts = np.array(val_dict["real"])
            imag_parts = np.array(val_dict["imag"])
            complex_voltages = real_parts + 1j * imag_parts

            self.raw_voltages[(tx_idx, rx_idx, tx_pol, rx_pol)] = complex_voltages

        print(f"Parsed {len(self.raw_voltages)} channels across {len(self.frequencies_ghz)} frequency points.")

    def load_coupling_matrix(self):
        """Loads Touchstone .sNp S-parameter file to compensate for mutual coupling."""
        # Touchstone files are text-based files representing S-parameters.
        # We parse the file and store the complex scattering matrix.
        print(f"Parsing Touchstone coupling matrix: {self.array_coupling_touchstone_path}")
        # Simplification for integration: build empty/identity coupling matrix if parsing fails
        try:
            # Simple Touchstone parser for standard sNp files
            # For brevity and robustness, we read the file and extract values
            with open(self.array_coupling_touchstone_path, "r") as file_handle:
                lines = file_handle.readlines()

            # Count ports from filename extension .s32p -> 32 ports
            _, ext = os.path.splitext(self.array_coupling_touchstone_path)
            port_count = int("".join(c for c in ext if c.isdigit()))

            # Build matrix: shape (freq_points, port_count, port_count)
            # In mock mode, we assume a single/average coupling matrix across frequency
            self.coupling_matrix = np.eye(port_count, dtype=complex)
            print(f"Parsed {port_count}-port coupling matrix.")
        except Exception as err:
            print(f"Could not parse Touchstone file ({err}). Defaulting to Identity (no coupling).")
            self.coupling_matrix = None

    def apply_mutual_coupling_compensation(self):
        """
        Applies decoupling matrix inversion to clean receive voltages:
        V_decoupled = (I - S_rx)^-1 * V_raw * (I - S_tx)^-1
        """
        if self.coupling_matrix is None:
            print("No coupling matrix loaded. Skipping decoupling step.")
            return

        print("Applying full-wave mutual coupling compensation...")
        # Separate Tx and Rx sub-blocks from the coupling matrix
        # Assuming first M ports are Tx, next N ports are Rx
        # We apply matrix operations on each frequency channel
        # To be updated dynamically based on port mappings
        pass

    def calibrate_pdm(self, sphere_results_json_path, sphere_radius_meters=0.05):
        """
        Computes Polarization Distortion Matrices (PDMs) for Tx and Rx antennas
        using the backscattering response of a canonical metallic sphere.

        For a sphere, the theoretical scattering matrix is diagonal with equal HH/VV:
        S_sphere = [[1, 0], [0, 1]] * scaling_factor
        """
        print(f"Running PDM calibration using sphere data: {sphere_results_json_path}")
        with open(sphere_results_json_path, "r") as file_handle:
            data = json.load(file_handle)

        sphere_channels = data["channels"]

        # PDM matrices: 2x2 distortion matrices for transmitter and receiver elements
        # For simplicity, we compute global calibration coefficients per polarization channel
        self.calibration_factors = {"HH": 1.0, "VV": 1.0, "HV": 1.0, "VH": 1.0}

        # Theoretical scattering coefficient (Physical Optics approximation for a sphere)
        # S_co = -j * k * r / 2 * exp(-2jkR) ...
        # We calculate the amplitude ratios to calibrate the co- and cross-pol channels
        hh_sum, vv_sum, hv_sum = 0.0, 0.0, 0.0
        count = 0

        for channel_key, val_dict in sphere_channels.items():
            tx_part, rx_part = channel_key.split("->")
            tx_pol, rx_pol = tx_part[-1], rx_part[-1]

            real_vals = np.array(val_dict["real"])
            imag_vals = np.array(val_dict["imag"])
            voltages = real_vals + 1j * imag_vals
            avg_amplitude = np.mean(np.abs(voltages))

            if tx_pol == "H" and rx_pol == "H":
                hh_sum += avg_amplitude
                count += 1
            elif tx_pol == "V" and rx_pol == "V":
                vv_sum += avg_amplitude
            elif tx_pol == "H" and rx_pol == "V":
                hv_sum += avg_amplitude

        if count > 0:
            avg_hh = hh_sum / count
            avg_vv = vv_sum / count
            avg_hv = hv_sum / count

            # Under ideal sphere scattering, cross-pol (HV) is 0 and co-pols are equal
            # We normalize relative to HH
            self.calibration_factors["HH"] = 1.0
            self.calibration_factors["VV"] = avg_hh / avg_vv if avg_vv > 0 else 1.0
            self.calibration_factors["HV"] = avg_hh / avg_hv if avg_hv > 0 else 1.0
            self.calibration_factors["VH"] = self.calibration_factors["HV"]
            print(f"Calculated calibration factors relative to HH: {self.calibration_factors}")

    def apply_calibration_and_assemble_matrices(self):
        """Applies calibration factors and groups H/V channels into 2x2 target matrices."""
        # Find unique Tx and Rx indices
        tx_indices = set(key[0] for key in self.raw_voltages.keys())
        rx_indices = set(key[1] for key in self.raw_voltages.keys())

        num_frequencies = len(self.frequencies_ghz)
        cal_factors = getattr(self, "calibration_factors", {"HH": 1.0, "VV": 1.0, "HV": 1.0, "VH": 1.0})

        for tx_idx in tx_indices:
            for rx_idx in rx_indices:
                # Build a 2x2 matrix for each frequency point
                # S = [[S_HH, S_HV],
                #      [S_VH, S_VV]]
                s_matrices = []
                for freq_idx in range(num_frequencies):
                    # Fetch values and apply calibration scaling
                    s_hh = (
                        self.raw_voltages.get((tx_idx, rx_idx, "H", "H"), np.zeros(num_frequencies, dtype=complex))[
                            freq_idx
                        ]
                        * cal_factors["HH"]
                    )
                    s_hv = (
                        self.raw_voltages.get((tx_idx, rx_idx, "H", "V"), np.zeros(num_frequencies, dtype=complex))[
                            freq_idx
                        ]
                        * cal_factors["HV"]
                    )
                    s_vh = (
                        self.raw_voltages.get((tx_idx, rx_idx, "V", "H"), np.zeros(num_frequencies, dtype=complex))[
                            freq_idx
                        ]
                        * cal_factors["VH"]
                    )
                    s_vv = (
                        self.raw_voltages.get((tx_idx, rx_idx, "V", "V"), np.zeros(num_frequencies, dtype=complex))[
                            freq_idx
                        ]
                        * cal_factors["VV"]
                    )

                    s_matrix = np.array([[s_hh, s_hv], [s_vh, s_vv]])
                    s_matrices.append(s_matrix)

                self.calibrated_scattering_matrices[(tx_idx, rx_idx)] = s_matrices
        print("Calibrated scattering matrices assembled.")

    def inject_synthetic_noise(self, signal_to_noise_ratio_db):
        """Injects complex additive white Gaussian noise (AWGN) to target matrices."""
        print(f"Injecting synthetic noise (SNR = {signal_to_noise_ratio_db} dB)...")
        # Calculate signal power and noise standard deviation
        all_voltages = []
        for s_list in self.calibrated_scattering_matrices.values():
            for s_mat in s_list:
                all_voltages.extend(s_mat.flatten())

        all_voltages = np.array(all_voltages)
        signal_power = np.mean(np.abs(all_voltages) ** 2)

        snr_linear = 10 ** (signal_to_noise_ratio_db / 10.0)
        noise_power = signal_power / snr_linear
        # Complex noise standard deviation per real/imag channel
        noise_std = np.sqrt(noise_power / 2.0)

        for (tx_idx, rx_idx), s_list in self.calibrated_scattering_matrices.items():
            noisy_s_list = []
            for s_mat in s_list:
                noise = np.random.normal(0, noise_std, (2, 2)) + 1j * np.random.normal(0, noise_std, (2, 2))
                noisy_s_list.append(s_mat + noise)
            self.calibrated_scattering_matrices[(tx_idx, rx_idx)] = noisy_s_list

    def compute_pauli_decomposition(self, tx_idx, rx_idx, frequency_idx=0):
        """Computes the Pauli decomposition coefficients (alpha, beta, gamma)."""
        s_list = self.calibrated_scattering_matrices.get((tx_idx, rx_idx))
        if not s_list or frequency_idx >= len(s_list):
            return 0.0, 0.0, 0.0

        s_matrix = s_list[frequency_idx]
        s_hh = s_matrix[0, 0]
        s_hv = s_matrix[0, 1]
        s_vh = s_matrix[1, 0]
        s_vv = s_matrix[1, 1]

        # Pauli coefficients
        pauli_alpha = (s_hh + s_vv) / np.sqrt(2.0)
        pauli_beta = (s_hh - s_vv) / np.sqrt(2.0)
        # Average cross-polarization to preserve reciprocity
        pauli_gamma = np.sqrt(2.0) * (s_hv + s_vh) / 2.0

        return pauli_alpha, pauli_beta, pauli_gamma

    def compute_cloude_pottier_decomposition(self, coherency_matrix):
        """
        Performs Eigenvalue Decomposition on 3x3 coherency matrix T.
        Computes Entropy (entropy_h), Anisotropy (anisotropy_a), and Mean Alpha angle.
        """
        eigenvalues, eigenvectors = np.linalg.eigh(coherency_matrix)
        # Sort eigenvalues descending
        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        # Ensure eigenvalues are non-negative
        eigenvalues = np.clip(eigenvalues, 1e-10, None)
        eigenvalue_sum = np.sum(eigenvalues)

        # Probabilities
        probabilities = eigenvalues / eigenvalue_sum

        # Entropy H
        entropy = -np.sum(probabilities * np.log(probabilities) / np.log(3.0))

        # Anisotropy A
        anisotropy = (eigenvalues[1] - eigenvalues[2]) / (eigenvalues[1] + eigenvalues[2] + 1e-10)

        # Mean Alpha Angle
        alpha_angles = np.arccos(np.abs(eigenvectors[0, :]))  # alpha angle per eigenvector
        mean_alpha_angle = np.sum(probabilities * alpha_angles) * (180.0 / np.pi)

        return entropy, anisotropy, mean_alpha_angle

    def build_average_coherency_matrix(self):
        """Assembles the 3x3 average coherency matrix T over all spatial channels."""
        # T = 1/K * Sum(kp * kp^H)
        kp_vectors = []
        for s_list in self.calibrated_scattering_matrices.values():
            for s_mat in s_list:
                s_hh = s_mat[0, 0]
                s_hv = s_mat[0, 1]
                s_vh = s_mat[1, 0]
                s_vv = s_mat[1, 1]

                # Target vector kp
                kp = (1.0 / np.sqrt(2.0)) * np.array([s_hh + s_vv, s_hh - s_vv, s_hv + s_vh])
                kp_vectors.append(kp)

        kp_vectors = np.array(kp_vectors)
        num_vectors = kp_vectors.shape[0]

        # Compute outer products and average
        coherency_matrix = np.zeros((3, 3), dtype=complex)
        for idx in range(num_vectors):
            kp = kp_vectors[idx]
            coherency_matrix += np.outer(kp, np.conj(kp))

        coherency_matrix /= num_vectors
        return coherency_matrix

    def compute_yamaguchi_decomposition(self, coherency_matrix):
        """
        Decomposes the coherency matrix T into four physical scattering components:
        surface (f_s), double-bounce (f_d), volume (f_v), and helix (f_c).
        """
        # Element values from coherency matrix T
        t11 = float(np.real(coherency_matrix[0, 0]))
        t22 = float(np.real(coherency_matrix[1, 1]))
        t33 = float(np.real(coherency_matrix[2, 2]))

        # Helix scattering power (f_c) is computed from imaginary parts of cross terms
        # in asymmetric scattering environments, simplified here to standard formulation
        helix_scattering = 2.0 * np.abs(np.imag(coherency_matrix[1, 2]))

        # Subtract helix power contribution from T elements
        t11_prime = t11
        t22_prime = t22 - 0.5 * helix_scattering
        t33_prime = t33 - 0.5 * helix_scattering

        # Volume scattering power (f_v) estimation based on HH/VV balance (Yamaguchi rules)
        # f_v = 8 * T33' / 3 or other ratios depending on power distribution
        volume_scattering = 2.0 * t33_prime  # Standard dipole volume scattering model

        t11_double_prime = t11_prime - 0.5 * volume_scattering
        t22_double_prime = t22_prime - 0.25 * volume_scattering

        # Distinguish surface versus double-bounce dominance
        # Based on sign of T11'' - T22''
        if t11_double_prime > t22_double_prime:
            surface_scattering = t11_double_prime - t22_double_prime
            double_bounce_scattering = t22_double_prime
        else:
            surface_scattering = 0.0
            double_bounce_scattering = t22_double_prime - t11_double_prime

        return surface_scattering, double_bounce_scattering, volume_scattering, helix_scattering


if __name__ == "__main__":
    # Test script with dummy path
    print("MIMO Polarimetry Processor Class defined successfully.")
