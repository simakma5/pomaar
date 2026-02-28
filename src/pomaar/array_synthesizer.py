from collections import defaultdict
import logging

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


def interleave_linear_arrays(tx_h: np.ndarray, tx_v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Swap every other element between two equal-length linear arrays."""
    tx_h = np.asarray(tx_h)
    tx_v = np.asarray(tx_v)
    if tx_h.shape != tx_v.shape:
        raise ValueError("tx_h and tx_v must have the same shape for interleaving")
    if tx_h.ndim == 1:
        tx_h = tx_h[:, None]
        tx_v = tx_v[:, None]

    out_h = tx_h.copy()
    out_v = tx_v.copy()
    if out_h.shape[0] > 1:
        idx = np.arange(out_h.shape[0]) % 2 == 1
        # swap using copy to avoid overwriting
        tmp = out_h[idx].copy()
        out_h[idx] = out_v[idx]
        out_v[idx] = tmp
    return out_h, out_v


class ArraySynthesizer:
    def __init__(self, name="2D polarimetric MIMO layout", highlight_overlaps=True):
        """Positions are in units of lambda/2, representing the electromagnetic phase centers of the antennas."""
        self.name = name
        self.highlight_overlaps = highlight_overlaps
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

        # Physical array
        self.tx_h = np.empty((0, 2))
        self.tx_v = np.empty((0, 2))
        self.rx_h = np.empty((0, 2))
        self.rx_v = np.empty((0, 2))
        self.tx_h_count = 0
        self.tx_v_count = 0
        self.tx_count = 0
        self.rx_count = 0

        # Virtual array
        self.v_hh = np.empty((0, 2))
        self.v_vv = np.empty((0, 2))
        self.v_hv = np.empty((0, 2))
        self.v_vh = np.empty((0, 2))
        self.virtual_element_count = 0

    # region public API
    def set_arrays(self, tx_h=[], rx_h=[], tx_v=[], rx_v=[], centre=False, centre_method="bounds"):
        """Sets the physical positions of the array elements (N x 2 arrays).

        The centre_method accepts two strategies:
        - "mean": Compute the arithmetic mean of all coordinates.
        - "bounds": Compute the mid-point of the bounding box (more consistent for non-uniform element distribution).
        """
        self.tx_h = np.array(tx_h) if len(tx_h) > 0 else np.empty((0, 2))
        self.tx_v = np.array(tx_v) if len(tx_v) > 0 else np.empty((0, 2))
        self.rx_h = np.array(rx_h) if len(rx_h) > 0 else np.empty((0, 2))
        self.rx_v = np.array(rx_v) if len(rx_v) > 0 else np.empty((0, 2))
        self.tx_h_count = self.tx_h.shape[0]
        self.tx_v_count = self.tx_v.shape[0]
        self.rx_h_count = self.rx_h.shape[0]
        self.rx_v_count = self.rx_v.shape[0]
        self.tx_count = self.tx_h_count + self.tx_v_count
        self.rx_count = self.rx_h_count + self.rx_v_count
        if centre:
            self._centre_physical(method=centre_method)
        self._compute_virtual()

    def plot_topology(self, xlim=None, ylim=None, marker_size=80):
        """Visualizes the physical and virtual arrays in 2D."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=False)
        for ax in [ax1, ax2]:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim)
        print_legend = True
        if any([not non_zero for non_zero in [self.tx_h_count, self.tx_v_count, self.rx_h_count, self.rx_v_count]]):
            self.highlight_overlaps = False
            print_legend = None
        # Physical layout
        if self.rx_h_count:
            ax1.scatter(
                self.rx_h[:, 0],
                self.rx_h[:, 1],
                s=marker_size,
                c="#0047AB",
                marker="s",
                label="Rx H" if print_legend else "Rx",
                edgecolors="k",
            )
        if self.rx_v_count:
            ax1.scatter(
                self.rx_v[:, 0],
                self.rx_v[:, 1],
                s=marker_size,
                c="#DC143C",
                marker="s",
                label="Rx V" if print_legend else "Rx",
                edgecolors="k",
            )
        if self.tx_h_count:
            ax1.scatter(
                self.tx_h[:, 0],
                self.tx_h[:, 1],
                s=marker_size,
                c="#0047AB",
                marker="^",
                label="Tx H" if print_legend else "Tx",
                edgecolors="k",
            )
        if self.tx_v_count:
            ax1.scatter(
                self.tx_v[:, 0],
                self.tx_v[:, 1],
                s=marker_size,
                c="#DC143C",
                marker="^",
                label="Tx V" if print_legend else "Tx",
                edgecolors="k",
            )
        ax1.set_title(
            f"Physical array ({self.tx_count} x {self.rx_count})",
            fontsize=14,
        )
        ax1.set_ylabel(r"Elevation ($\lambda/2$)", fontsize=12)
        ax1.legend(loc="best")
        ax1.grid(True, linestyle="--", alpha=0.5)
        # Virtual array
        offset = 0  # optional offset to separate overlapping channels
        if self.v_hh.size:
            ax2.scatter(
                self.v_hh[:, 0],
                self.v_hh[:, 1] - offset,
                s=marker_size,
                c="#0047AB",
                marker="+",
                label="HH" if print_legend else None,
            )
        if self.v_vv.size:
            ax2.scatter(
                self.v_vv[:, 0],
                self.v_vv[:, 1] + offset,
                s=marker_size,
                c="#DC143C",
                marker="+",
                label="VV" if print_legend else None,
            )
        if self.v_hv.size:
            ax2.scatter(
                self.v_hv[:, 0] - offset,
                self.v_hv[:, 1] - offset,
                s=marker_size,
                c="#2E8B57",
                marker="x",
                label="HV" if print_legend else None,
            )
        if self.v_vh.size:
            ax2.scatter(
                self.v_vh[:, 0] + offset,
                self.v_vh[:, 1] + offset,
                s=marker_size,
                c="#8A2BE2",
                marker="x",
                label="VH" if print_legend else None,
            )
        # Highlight calibration overlaps
        calibration_overlaps, redundant_overlaps = self._analyze_calibration_overlaps()
        total_virtual_elements = (
            self.tx_count * self.rx_count - calibration_overlaps.shape[0] - redundant_overlaps.shape[0]
        )
        if self.highlight_overlaps:
            if calibration_overlaps.size:
                ax2.scatter(
                    calibration_overlaps[:, 0],
                    calibration_overlaps[:, 1],
                    s=1.5 * marker_size,
                    facecolors="none",
                    edgecolors="gold",
                    linewidth=1.5,
                    label="Calibration overlap",
                )
                overlap_str = self._format_coords(calibration_overlaps)
                self.logger.info(f"{len(calibration_overlaps)} calibration overlaps found: {overlap_str}")

                if redundant_overlaps.size:
                    ax2.scatter(
                        redundant_overlaps[:, 0],
                        redundant_overlaps[:, 1],
                        s=1.5 * marker_size,
                        facecolors="none",
                        edgecolors="red",
                        linewidth=1.5,
                        label="Redundant overlap",
                    )
                    redundant_str = self._format_coords(redundant_overlaps)
                    self.logger.warning(f"{len(redundant_overlaps)} redundant overlaps found: {redundant_str}")
            else:
                print("No overlaps found.")

        ax2.set_title(f"Virtual array ({total_virtual_elements} unique elements)", fontsize=14)
        ax2.set_xlabel(r"Azimuth ($\lambda/2$)", fontsize=12)
        ax2.set_ylabel(r"Elevation ($\lambda/2$)", fontsize=12)
        if print_legend:
            ax2.legend(loc="best", ncol=1)
        ax2.grid(True, linestyle="--", alpha=0.5)
        fig.suptitle(self.name, fontsize=16)
        plt.tight_layout()
        plt.show()

    # endregion
    # region private methods

    def _analyze_calibration_overlaps(self):
        """
        Identifies channel overlaps at each virtual position, distinguishing between
        calibration and redundant ones. Utilizes a binning approach to group nearby
        virtual elements, accounting for minor numerical discrepancies.
        """
        resolution = 1e-2  # Main tuning knob for adjusting the binning sensitivity
        occupied_bins = defaultdict(set)
        for channel, positions in [("hh", self.v_hh), ("vv", self.v_vv), ("hv", self.v_hv), ("vh", self.v_vh)]:
            if not positions.size:
                continue
            binned_positions = np.floor(positions / resolution).astype(int)
            for position in binned_positions:
                occupied_bins[tuple(position)].add(channel)
        calibration_overlaps, redundant_overlaps = [], []
        for position, channels in occupied_bins.items():
            if len(channels) < 2:
                continue
            real_coord = [coord * resolution for coord in position]
            if any(pair.issubset(channels) for pair in [{"hh", "vv"}, {"hv", "vh"}]):
                calibration_overlaps.append(real_coord)
            else:
                redundant_overlaps.append(real_coord)

        return np.array(sorted(calibration_overlaps)), np.array(sorted(redundant_overlaps))

    def _centre_physical(self, method: str = "bounds") -> np.ndarray:
        """Translate every physical coordinate so that the layout is centred."""
        all_positions = np.vstack([arr for arr in (self.tx_h, self.tx_v, self.rx_h, self.rx_v) if arr.size])
        if all_positions.size == 0:
            centre = np.zeros(2, dtype=float)
        if method.lower() == "mean":
            centre = all_positions.mean(axis=0)
        elif method.lower() == "bounds":
            minimum = all_positions.min(axis=0)
            maximum = all_positions.max(axis=0)
            centre = (minimum + maximum) / 2.0
        else:
            raise ValueError(f"Unknown centring method '{method}'")
        for arr_name in ("tx_h", "tx_v", "rx_h", "rx_v"):
            arr = getattr(self, arr_name)
            if arr.size:
                setattr(self, arr_name, arr - centre)
        return centre

    def _compute_virtual(self):
        """Computes the spatial convolution (MIMO virtual array) in 2D."""

        def convolve(tx: np.ndarray, rx: np.ndarray):
            """
            Computes virtual positions via the sum of transmit and receive vectors.

            This uses the 'sum convention' (x_v = x_t + x_r) rather than the physical
            midpoint. While the phase centres physically lie at (x_t + x_r) / 2,
            omitting the division by 2 allows the virtual array to be treated as
            a receive-only aperture receiving a one-way plane wave, maintaining
            the standard d = lambda/2 Nyquist spacing for the virtual elements.

            (N_tx, 1, 2) + (1, N_rx, 2) -> (N_tx, N_rx, 2) -> (N_tx * N_rx, 2)
            """
            return (tx[:, None, :] + rx[None, :, :]).reshape(-1, 2) if tx.size and rx.size else np.empty((0, 2))

        self.v_hh = convolve(self.tx_h, self.rx_h)
        self.v_vv = convolve(self.tx_v, self.rx_v)
        self.v_hv = convolve(self.tx_h, self.rx_v)
        self.v_vh = convolve(self.tx_v, self.rx_h)

    def _format_coords(self, coords: np.ndarray, max_items=None):
        """Internal helper to prettify coordinate lists for logging."""
        if not coords.size:
            return "[]"
        # 'f' for fixed-point notation; 'g' for scientific notation
        formatted_list = [f"({x:.1f}, {y:.1f})" for x, y in coords]
        if max_items is not None and len(formatted_list) > max_items:
            return f"{', '.join(formatted_list[:max_items])} ... {formatted_list[-1]}"
        return ", ".join(formatted_list)

    # endregion
    # region example


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # Receiver array (ULA)
    rx_h = []
    rx_v = []
    # Interleaved pattern: H, V, H, V...
    n_rx_chips = 4
    elements_per_chip = 4
    for x in range(n_rx_chips * elements_per_chip):
        pos = [x, 0]
        if x % 2 == 0:
            rx_h.append(pos)
        else:
            rx_v.append(pos)

    # Transmitter array (mild non-uniformity in elevation)
    # Azimuth Tx
    tx_h_az = [[-4, 0], [20, 0]]
    tx_v_az = [[-3, 0], [19, 0]]  # Slight offset to create overlap with Rx shift
    # Elevation Tx
    tx_h_el = [[-4, 5]]
    tx_v_el = [[20, 5]]

    tx_h = np.array(tx_h_az + tx_h_el)
    tx_v = np.array(tx_v_az + tx_v_el)

    # Initialize and run
    radar = ArraySynthesizer(name="test layout")
    radar.set_arrays(tx_h, tx_v, rx_h, rx_v)
    radar.plot_topology()
