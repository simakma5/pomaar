#!/usr/bin/env python3
import argparse

from _layout_common import (
    add_unit_cell_arguments,
    layout_filename,
    unit_cell_metadata,
    write_layout,
)


def main():
    parser = argparse.ArgumentParser(description="Generate a Polarimetric ULA MIMO Layout YAML.")
    parser.add_argument(
        "--rx-count", type=int, default=4, help="Number of Rx elements (default: 4)"
    )
    parser.add_argument(
        "--tx-count", type=int, default=4, help="Number of Tx elements (default: 4)"
    )
    parser.add_argument(
        "--centre-freq",
        "--center-freq",
        type=float,
        default=79.0,
        help=(
            "Design frequency in GHz for the default lambda-based spacings (default: 79.0); "
            "the builder takes the simulation frequencies from the unit cell"
        ),
    )
    parser.add_argument(
        "--rx-spacing",
        type=float,
        default=None,
        help="Spacing between Rx elements in mm (default: 0.5 * lambda_0)",
    )
    parser.add_argument(
        "--tx-spacing",
        type=float,
        default=None,
        help="Spacing between Tx elements in mm (default: rx_count * rx_spacing)",
    )
    parser.add_argument(
        "--tx-offset-y",
        type=float,
        default=12.0,
        help="Total Y separation between Rx and Tx rows in mm (default: 12.0)",
    )
    parser.add_argument(
        "--pcb-margin", type=float, default=0.0, help="Additional PCB margin in mm (default: 0.0)"
    )

    add_unit_cell_arguments(parser)
    args = parser.parse_args()

    # Dynamically calculate spacing based on the centre frequency if not specified
    lambda_0 = 299.792458 / args.centre_freq
    if args.rx_spacing is None:
        args.rx_spacing = round(0.5 * lambda_0, 2)

    if args.tx_spacing is None:
        args.tx_spacing = args.rx_count * args.rx_spacing

    # Element positions are centred around the array centroid.
    # Rx row sits at y = -txOffsetY/2, Tx row at y = +txOffsetY/2.
    rx_y = -args.tx_offset_y / 2.0
    tx_y = args.tx_offset_y / 2.0

    # Every Rx is one dual-polarized element (polarization: hv) built from a dual-port unit cell;
    # two single-pol copies at one position would overlap. The Tx elements are single-pol copies
    # of the source design, assumed V-polarized with its feed towards -y at yaw 0, and H is made by
    # turning it. Their yaws keep the pair sign s_H * s_V = +1 (V yaw 180 -> -y, H yaw 90 -> -x),
    # so no cross-pol channel is measured with a 180 deg flip relative to the co-pol ones.
    elements = []
    rx_offset_x = (args.rx_count - 1) * args.rx_spacing / 2.0
    for i in range(args.rx_count):
        x_pos = i * args.rx_spacing - rx_offset_x
        rx_idx_offset = f"({i} - ({args.rx_count} - 1) / 2.0)" if args.rx_count > 1 else "0"
        elements.append(
            {
                "label": f"Rx_{i + 1}",
                "role": "Rx",
                "position": [round(x_pos, 4), round(rx_y, 4), 0.0],
                "position_expression": [f"{rx_idx_offset} * rxSpacing", "-txOffsetY / 2"],
                "polarization": "hv",
                "yaw": 0.0,
            }
        )

    # Alternating polarizations for Transmitters
    tx_offset_x = (args.tx_count - 1) * args.tx_spacing / 2.0
    for i in range(args.tx_count):
        x_pos = i * args.tx_spacing - tx_offset_x
        tx_idx_offset = f"({i} - ({args.tx_count} - 1) / 2.0)" if args.tx_count > 1 else "0"
        is_v_pol = i % 2 == 0
        elements.append(
            {
                "label": f"Tx_{i + 1}",
                "role": "Tx",
                "position": [round(x_pos, 4), round(tx_y, 4), 0.0],
                "position_expression": [f"{tx_idx_offset} * txSpacing", "txOffsetY / 2"],
                "polarization": "v" if is_v_pol else "h",
                "yaw": 180.0 if is_v_pol else 90.0,
            }
        )

    layout_data = {
        "metadata": {
            **unit_cell_metadata(parser, args),
            "variables": {
                "rxSpacing": f"{args.rx_spacing:.4f}mm",
                "txSpacing": f"{args.tx_spacing:.4f}mm",
                "txOffsetY": f"{args.tx_offset_y:.4f}mm",
                "pcbMargin": f"{args.pcb_margin:.4f}mm",
            },
            "board": {
                "width_formula": f"max(({args.rx_count} - 1) * rxSpacing, ({args.tx_count} - 1) * txSpacing) + 2 * (unitCellExtent + pcbMargin)",
                "length_formula": "txOffsetY + 2 * (unitCellExtent + pcbMargin)",
            },
        },
        "elements": elements,
    }
    if "hv" not in layout_data["metadata"].get("unit_cells", {}):
        print(
            "Note: the dual-polarized Rx elements need a dual-port unit cell; regenerate with "
            "--extra-unit-cell hv=<Design> before building."
        )

    output_path = layout_filename(args, "polarimetric_ula")

    print(
        f"Generated layout with {args.rx_count} Rx (dual-pol) and {args.tx_count} Tx (alternating V/H) elements."
    )
    print(f"  Centre frequency: {args.centre_freq} GHz (lambda_0 = {lambda_0:.2f} mm)")
    print(f"  Rx Spacing (Dual-Pol): {args.rx_spacing} mm")
    print(f"  Tx Spacing (Alternating V/H): {args.tx_spacing:.2f} mm")
    print(f"  Rx row at y = {rx_y:.2f} mm, Tx row at y = {tx_y:.2f} mm (centred)")
    write_layout(layout_data, output_path, __file__)


if __name__ == "__main__":
    main()
