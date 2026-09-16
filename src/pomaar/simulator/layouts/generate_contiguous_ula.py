#!/usr/bin/env python3
import argparse
from pathlib import Path
import yaml


def main():
    parser = argparse.ArgumentParser(description="Generate a Contiguous ULA MIMO Layout YAML.")
    parser.add_argument("--rx-count", type=int, default=4, help="Number of Rx elements (default: 4)")
    parser.add_argument("--tx-count", type=int, default=4, help="Number of Tx elements (default: 4)")
    parser.add_argument(
        "--centre-freq", "--center-freq", type=float, default=79.0, help="Centre frequency in GHz (default: 79.0)"
    )
    parser.add_argument("--bandwidth", type=float, default=4.0, help="Sweep bandwidth in GHz (default: 4.0)")
    parser.add_argument(
        "--rx-spacing", type=float, default=None, help="Spacing between Rx elements in mm (default: 0.5 * lambda_0)"
    )
    parser.add_argument(
        "--tx-spacing",
        type=float,
        default=None,
        help="Spacing between Tx elements in mm (default: rx_count * rx_spacing)",
    )
    parser.add_argument(
        "--tx-offset-y", type=float, default=10.0, help="Total Y separation between Rx and Tx rows in mm (default: 10.0)"
    )
    parser.add_argument(
        "--pcb-margin", type=float, default=0.0, help="Additional PCB margin in mm (default: 0.0)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="contiguous_ula.yaml",
        help="Output YAML filename (appended with _<tx_count>x<rx_count>, default: contiguous_ula.yaml)",
    )

    args = parser.parse_args()

    # Dynamically calculate spacing based on centre frequency if not specified
    lambda_0 = 299.792458 / args.centre_freq
    if args.rx_spacing is None:
        args.rx_spacing = round(0.5 * lambda_0, 2)

    if args.tx_spacing is None:
        args.tx_spacing = args.rx_count * args.rx_spacing

    # Element positions are centered around the array centroid.
    # Rx row sits at y = -txOffsetY/2, Tx row at y = +txOffsetY/2.
    rx_y = -args.tx_offset_y / 2.0
    tx_y = args.tx_offset_y / 2.0

    # Construct Rx Elements (dense ULA elements)
    elements = []
    rx_offset_x = (args.rx_count - 1) * args.rx_spacing / 2.0
    for i in range(args.rx_count):
        x_pos = i * args.rx_spacing - rx_offset_x
        rx_idx_offset = f"({i} - ({args.rx_count} - 1) / 2.0)" if args.rx_count > 1 else "0"
        elements.append(
            {
                "label": f"Rx_{i+1}",
                "role": "Rx",
                "position": [round(x_pos, 4), round(rx_y, 4), 0.0],
                "position_expression": [f"{rx_idx_offset} * rxSpacing", "-txOffsetY / 2", "0mm"],
                "polarization": "v",
                "yaw": 0.0,
            }
        )

    # Construct Tx Elements (sparse ULA elements)
    tx_offset_x = (args.tx_count - 1) * args.tx_spacing / 2.0
    for i in range(args.tx_count):
        x_pos = i * args.tx_spacing - tx_offset_x
        tx_idx_offset = f"({i} - ({args.tx_count} - 1) / 2.0)" if args.tx_count > 1 else "0"
        elements.append(
            {
                "label": f"Tx_{i+1}",
                "role": "Tx",
                "position": [round(x_pos, 4), round(tx_y, 4), 0.0],
                "position_expression": [f"{tx_idx_offset} * txSpacing", "txOffsetY / 2", "0mm"],
                "polarization": "v",
                "yaw": 180.0,
            }
        )

    layout_data = {
        "metadata": {
            "topology": "contiguous_ula",
            "center_frequency_ghz": args.centre_freq,
            "bandwidth_ghz": args.bandwidth,
            "variables": {
                "rxSpacing": f"{args.rx_spacing:.4f}mm",
                "txSpacing": f"{args.tx_spacing:.4f}mm",
                "txOffsetY": f"{args.tx_offset_y:.4f}mm",
                "pcbMargin": f"{args.pcb_margin:.4f}mm",
            },
            "board": {
                "width_formula": f"max(({args.rx_count} - 1) * rxSpacing, ({args.tx_count} - 1) * txSpacing) + 2 * (unitCellExtentX + pcbMargin)",
                "length_formula": "txOffsetY + 2 * (unitCellExtentY + pcbMargin)",
            },
        },
        "elements": elements,
    }

    output_path = Path(args.output)
    tag = f"_{args.tx_count}x{args.rx_count}"
    if not output_path.stem.endswith(tag):
        output_path = output_path.with_name(f"{output_path.stem}{tag}{output_path.suffix}")

    if not output_path.is_absolute():
        output_path = Path(__file__).resolve().parent / output_path

    with open(output_path, "w", encoding="utf-8") as f:
        yaml.dump(layout_data, f, sort_keys=False)

    print(f"Generated layout with {args.rx_count} Rx and {args.tx_count} Tx elements.")
    print(f"  Centre Frequency: {args.centre_freq} GHz (lambda_0 = {lambda_0:.2f} mm)")
    print(f"  Rx Spacing: {args.rx_spacing} mm")
    print(f"  Tx Spacing: {args.tx_spacing:.2f} mm")
    print(f"  Rx row at y = {rx_y:.2f} mm, Tx row at y = {tx_y:.2f} mm (centered)")
    print(f"Layout written to: {output_path}")
    print(f"\nTo synthesize in HFSS, run:")
    print(
        f"  hfss_array_builder <project_path> <source_design_name> {output_path} --centre-freq {args.centre_freq} --bandwidth {args.bandwidth}"
    )


if __name__ == "__main__":
    main()
