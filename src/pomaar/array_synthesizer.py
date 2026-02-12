from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

# TODO: Visualize topologies from the literature (Tinti, TI manuals, etc.) to verify the code and provide examples.


class PoarimetricArraySynthesizer:
    def __init__(self, name="2D polarimetric MIMO layout"):
        self.name = name
        # Positions are in units of lambda/2
        self.tx_h = np.empty((0, 2))
        self.tx_v = np.empty((0, 2))
        self.rx_h = np.empty((0, 2))
        self.rx_v = np.empty((0, 2))

        # Virtual arrays
        self.v_hh = np.empty((0, 2))
        self.v_vv = np.empty((0, 2))
        self.v_hv = np.empty((0, 2))
        self.v_vh = np.empty((0, 2))

    def set_arrays(self, tx_h: list, tx_v: list, rx_h: list, rx_v: list):
        """Sets the physical positions of the array elements (N x 2 arrays)."""
        self.tx_h = np.array(tx_h) if len(tx_h) > 0 else np.empty((0, 2))
        self.tx_v = np.array(tx_v) if len(tx_v) > 0 else np.empty((0, 2))
        self.rx_h = np.array(rx_h) if len(rx_h) > 0 else np.empty((0, 2))
        self.rx_v = np.array(rx_v) if len(rx_v) > 0 else np.empty((0, 2))
        self._compute_virtual()

    def _compute_virtual(self):
        """Computes the spatial convolution (MIMO virtual array) in 2D."""

        def convolve(tx: np.ndarray, rx: np.ndarray):
            """(N_tx, 1, 2) + (1, N_rx, 2) -> (N_tx, N_rx, 2) -> (N_tx * N_rx, 2)"""
            return (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2) if tx.size and rx.size else np.empty((0, 2))

        self.v_hh = convolve(self.tx_h, self.rx_h)
        self.v_vv = convolve(self.tx_v, self.rx_v)
        self.v_hv = convolve(self.tx_h, self.rx_v)
        self.v_vh = convolve(self.tx_v, self.rx_h)

    def analyze_calibration_overlaps(self):
        """
        Identifies positions where any two virtual channels coincide, distinguishing between intended calibration
        overlaps (co-polar and cross-polar pairs) and redundant ones.
        """
        # Determine channels at each binned virtual position
        resolution = 1e-2
        occupied_bins = defaultdict(set)
        for channel, positions in [("hh", self.v_hh), ("vv", self.v_vv), ("hv", self.v_hv), ("vh", self.v_vh)]:
            if not positions.size:
                continue
            binned_positions = np.floor(positions / resolution).astype(int)
            for position in binned_positions:
                occupied_bins[tuple(position)].add(channel)
        # Distinguish between useful and redundant overlaps
        calibration_overlaps, redundant_overlaps = [], []
        for position, channels in occupied_bins.items():
            if len(channels) < 2:
                continue
            real_coord = [coord * resolution for coord in position]
            if any(pair.issubset(channels) for pair in [{"hh", "vv"}, {"hv", "vh"}]):
                calibration_overlaps.append(real_coord)
            else:
                redundant_overlaps.append(real_coord)

        return sorted(calibration_overlaps), sorted(redundant_overlaps)

    def plot_topology(self):
        """Visualizes the physical and virtual arrays in 2D."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

        # --- 1. Physical layout ---
        # Plot Rx
        if self.rx_h.size:
            ax1.scatter(self.rx_h[:, 0], self.rx_h[:, 1], s=80, c="#0047AB", marker="s", label="Rx H", edgecolors="k")
        if self.rx_v.size:
            ax1.scatter(self.rx_v[:, 0], self.rx_v[:, 1], s=80, c="#DC143C", marker="s", label="Rx V", edgecolors="k")

        # Plot Tx
        if self.tx_h.size:
            ax1.scatter(self.tx_h[:, 0], self.tx_h[:, 1], s=80, c="#0047AB", marker="^", label="Tx H", edgecolors="k")
        if self.tx_v.size:
            ax1.scatter(self.tx_v[:, 0], self.tx_v[:, 1], s=80, c="#DC143C", marker="^", label="Tx V", edgecolors="k")

        ax1_title = f"{self.name}\nPhysical array" if self.name else "Physical array"
        ax1.set_title(ax1_title, fontsize=14)
        ax1.set_ylabel(r"Elevation ($\lambda/2$)", fontsize=12)
        ax1.legend(loc="center right")
        ax1.grid(True, linestyle="--", alpha=0.5)

        # --- 2. Virtual array layout ---
        # Co-polar
        if self.v_hh.size:
            ax2.scatter(self.v_hh[:, 0], self.v_hh[:, 1], s=80, c="#0047AB", marker="+", label="HH")
        if self.v_vv.size:
            ax2.scatter(self.v_vv[:, 0], self.v_vv[:, 1], s=80, c="#DC143C", marker="+", label="VV")

        # Cross-polar (with the option of a slight offset for visibility if perfectly overlapping)
        offset = 0  # 0.05
        if self.v_hv.size:
            ax2.scatter(self.v_hv[:, 0] - offset, self.v_hv[:, 1] - offset, s=80, c="#2E8B57", marker="x", label="HV")
        if self.v_vh.size:
            ax2.scatter(self.v_vh[:, 0] + offset, self.v_vh[:, 1] + offset, s=80, c="#8A2BE2", marker="x", label="VH")

        # Highlight calibration overlaps
        calibration_overlaps, redundant_overlaps = self.analyze_calibration_overlaps()
        if calibration_overlaps:
            overlap_arr = np.array(list(calibration_overlaps))
            ax2.scatter(
                overlap_arr[:, 0],
                overlap_arr[:, 1],
                s=200,
                facecolors="none",
                edgecolors="gold",
                linewidth=1.5,
                label="Calibration overlap",
            )
            print(f"Calibration overlaps found at virtual positions: {calibration_overlaps}")  # TODO: log.info

            if redundant_overlaps:
                redundant_arr = np.array(list(redundant_overlaps))
                ax2.scatter(
                    redundant_arr[:, 0],
                    redundant_arr[:, 1],
                    s=200,
                    facecolors="none",
                    edgecolors="red",
                    linewidth=1.5,
                    label="Redundant overlap",
                )
                print(f"Redundant overlaps found at virtual positions: {redundant_overlaps}")  # TODO: log.warn
        else:
            print("No overlaps found.")

        ax2.set_title("Virtual array aperture", fontsize=14)
        ax2.set_xlabel(r"Azimuth ($\lambda/2$)", fontsize=12)
        ax2.set_ylabel(r"Elevation ($\lambda/2$)", fontsize=12)
        ax2.legend(loc="center right", ncol=1)
        ax2.grid(True, linestyle="--", alpha=0.5)

        plt.tight_layout()
        plt.show()


# Warning: Physical element positions represent antennas' electromagnetic phase centres -- not geometrical!
if __name__ == "__main__":
    # --- Example: TI cascade-style layout (2D) ---
    # [x, y] in lambda/2 units.

    # 1. Receiver array (ULA)
    rx_h_pos = []
    rx_v_pos = []
    # Interleaved pattern: H, V, H, V...
    n_rx_chips = 4
    elements_per_chip = 4
    for x in range(n_rx_chips * elements_per_chip):
        pos = [x, 0]
        if x % 2 == 0:
            rx_h_pos.append(pos)
        else:
            rx_v_pos.append(pos)

    # 2. Transmitter array (Sparse 2D)
    # To get elevation resolution, we need Tx elements at different y-coordinates.
    # To get azimuth resolution, we need Tx elements spaced widely in x.
    # Azimuth Tx (at y = 0)
    tx_h_az = [[-4, 0], [20, 0]]
    tx_v_az = [[-3, 0], [19, 0]]  # Slight offset to create overlap with Rx shift
    # Elevation Tx (at y != 0)
    tx_h_el = [[-4, 5]]
    tx_v_el = [[20, 5]]

    # Combine lists
    tx_h_pos = np.array(tx_h_az + tx_h_el)
    tx_v_pos = np.array(tx_v_az + tx_v_el)

    # Initialize and run
    radar = PoarimetricArraySynthesizer(name="test layout")
    radar.set_arrays(tx_h_pos, tx_v_pos, rx_h_pos, rx_v_pos)
    radar.plot_topology()
