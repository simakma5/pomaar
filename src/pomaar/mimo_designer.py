import matplotlib.pyplot as plt
import numpy as np


class PolarimetricMIMO:
    def __init__(self, name="MIMO Design"):
        self.name = name
        # Positions are in units of lambda/2
        self.tx_h = []
        self.tx_v = []
        self.rx_h = []
        self.rx_v = []

        # Virtual arrays
        self.v_hh = []
        self.v_vv = []
        self.v_hv = []
        self.v_vh = []

    def set_arrays(self, tx_h, tx_v, rx_h, rx_v):
        """Sets the physical positions of the array elements."""
        self.tx_h = np.array(tx_h)
        self.tx_v = np.array(tx_v)
        self.rx_h = np.array(rx_h)
        self.rx_v = np.array(rx_v)
        self._compute_virtual()

    def _compute_virtual(self):
        """Computes the spatial convolution (MIMO virtual array)."""
        # Virtual = Tx + Rx
        # We use standard broadcasting to get all combinations
        self.v_hh = (self.tx_h[:, None] + self.rx_h).flatten() if len(self.tx_h) and len(self.rx_h) else np.array([])
        self.v_vv = (self.tx_v[:, None] + self.rx_v).flatten() if len(self.tx_v) and len(self.rx_v) else np.array([])

        # Cross-polar terms
        self.v_hv = (self.tx_h[:, None] + self.rx_v).flatten() if len(self.tx_h) and len(self.rx_v) else np.array([])
        self.v_vh = (self.tx_v[:, None] + self.rx_h).flatten() if len(self.tx_v) and len(self.rx_h) else np.array([])

    def analyze_calibration_overlaps(self):
        """Identifies virtual positions where HV and VH coincide."""
        # Find intersection between HV and VH sets
        # Rounding is necessary to avoid floating point equality issues
        hv_set = set(np.round(self.v_hv, 2))
        vh_set = set(np.round(self.v_vh, 2))

        overlaps = hv_set.intersection(vh_set)
        return sorted(list(overlaps))

    def plot_topology(self):
        """Visualizes the physical and virtual arrays."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # --- 1. Physical Layout ---
        # Plot Rx (H and V)
        ax1.scatter(self.rx_h, np.zeros_like(self.rx_h), c="blue", marker="s", label="Rx H", s=100, edgecolors="k")
        ax1.scatter(self.rx_v, np.zeros_like(self.rx_v), c="red", marker="s", label="Rx V", s=100, edgecolors="k")

        # Plot Tx (H and V) - Offset slightly in Y for visibility
        ax1.scatter(self.tx_h, np.ones_like(self.tx_h) * 0.5, c="cyan", marker="^", label="Tx H", s=120, edgecolors="k")
        ax1.scatter(
            self.tx_v, np.ones_like(self.tx_v) * 0.5, c="magenta", marker="^", label="Tx V", s=120, edgecolors="k"
        )

        ax1.set_title(f"Physical Array Layout: {self.name}", fontsize=14)
        ax1.set_yticks([0, 0.5])
        ax1.set_yticklabels(["Rx Line", "Tx Line"])
        ax1.legend(loc="upper right")
        ax1.grid(True, linestyle="--", alpha=0.5)

        # --- 2. Virtual Array Layout ---
        # We plot these with some transparency to show density/overlap
        y_offset = 0

        # Co-polar
        ax2.scatter(self.v_hh, np.zeros_like(self.v_hh) + 0.2, c="blue", alpha=0.6, label="Virtual HH", s=60)
        ax2.scatter(self.v_vv, np.zeros_like(self.v_vv) - 0.2, c="red", alpha=0.6, label="Virtual VV", s=60)

        # Cross-polar
        ax2.scatter(
            self.v_hv, np.zeros_like(self.v_hv) + 0.05, c="purple", marker="x", alpha=0.8, label="Virtual HV", s=80
        )
        ax2.scatter(
            self.v_vh, np.zeros_like(self.v_vh) - 0.05, c="green", marker="+", alpha=0.8, label="Virtual VH", s=100
        )

        # Highlight Calibration Overlaps
        overlaps = self.analyze_calibration_overlaps()
        if overlaps:
            ax2.scatter(
                overlaps,
                np.zeros_like(overlaps),
                s=300,
                facecolors="none",
                edgecolors="gold",
                linewidth=2,
                label="Calib. Overlap (HV=VH)",
            )
            print(f"Calibration Overlaps found at virtual indices: {overlaps}")
        else:
            print("No HV/VH calibration overlaps found.")

        ax2.set_title("Virtual Array Aperture", fontsize=14)
        ax2.set_xlabel(r"Position ($\lambda/2$)", fontsize=12)
        ax2.set_yticks([])
        ax2.legend(loc="upper right", ncol=2)
        ax2.grid(True, axis="x", linestyle="--", alpha=0.5)

        plt.tight_layout()
        plt.show()


# --- Configuration Section: Interleaved Design Example ---
if __name__ == "__main__":
    # Let's assume a "unit" of 4 Rx elements (2H, 2V) interleaved
    # Rx Pattern: H V H V -> positions 0, 1, 2, 3
    rx_h_unit = [0, 2]
    rx_v_unit = [1, 3]

    # Create a larger array by repeating this unit (e.g., 2 chips)
    # Chip 1 at 0, Chip 2 at 4 (seamless abutment)
    rx_h_pos = np.concatenate([rx_h_unit, np.array(rx_h_unit) + 4])
    rx_v_pos = np.concatenate([rx_v_unit, np.array(rx_v_unit) + 4])

    # Tx Design:
    # To get overlaps, we need specific spacing.
    # If Rx array is length 8 (0 to 7), standard MIMO Tx spacing might be 8.
    # To force HV/VH overlap, we can use the symmetry trick.
    # Let's place H and V transmitters symmetrically around the center.
    tx_h_pos = [-4, 12]  # Wide aperture for resolution
    tx_v_pos = [-3, 11]  # Shifted by 1 unit relative to Tx_H (since Rx H/V are shifted by 1)

    # Initialize and run
    radar = PolarimetricMIMO(name="Interleaved 2-Chip Prototype")
    radar.set_arrays(tx_h_pos, tx_v_pos, rx_h_pos, rx_v_pos)
    radar.plot_topology()
