#!/usr/bin/env python3
import argparse
import sys
import yaml

def main():
    parser = argparse.ArgumentParser(description="Generate a Contiguous ULA MIMO Layout YAML.")
    parser.add_argument("--rx-count", type=int, default=4, help="Number of Rx elements (default: 4)")
    parser.add_argument("--tx-count", type=int, default=4, help="Number of Tx elements (default: 4)")
    parser.add_argument("--centre-freq", "--center-freq", type=float, default=79.0, help="Centre frequency in GHz (default: 79.0)")
    parser.add_argument("--bandwidth", type=float, default=4.0, help="Sweep bandwidth in GHz (default: 4.0)")
    parser.add_argument("--rx-spacing", type=float, default=None, help="Spacing between Rx elements in mm (default: 0.5 * lambda_0)")
    parser.add_argument("--tx-spacing", type=float, default=None, help="Spacing between Tx elements in mm (default: rx_count * rx_spacing)")
    parser.add_argument("--tx-offset-y", type=float, default=10.0, help="Y-offset of the Tx array in mm (default: 10.0)")
    parser.add_argument("-o", "--output", default="contiguous_ula_layout.yaml", help="Output YAML filename (default: contiguous_ula_layout.yaml)")

    args = parser.parse_args()

    # Dynamically calculate spacing based on centre frequency if not specified
    lambda_0 = 299.792458 / args.centre_freq
    if args.rx_spacing is None:
        args.rx_spacing = round(0.5 * lambda_0, 2)

    if args.tx_spacing is None:
        args.tx_spacing = args.rx_count * args.rx_spacing

    # Construct Rx Elements (dense ULA elements)
    elements = []
    rx_offset_x = (args.rx_count - 1) * args.rx_spacing / 2.0
    for i in range(args.rx_count):
        x_pos = i * args.rx_spacing - rx_offset_x
        elements.append({
            "label": f"Rx_{i+1}",
            "role": "Rx",
            "pos": [x_pos, 0.0, 0.0],
            "polarization": "v",
            "yaw": 0.0
        })

    # Construct Tx Elements (sparse ULA elements)
    tx_offset_x = (args.tx_count - 1) * args.tx_spacing / 2.0
    for i in range(args.tx_count):
        x_pos = i * args.tx_spacing - tx_offset_x
        elements.append({
            "label": f"Tx_{i+1}",
            "role": "Tx",
            "pos": [x_pos, args.tx_offset_y, 0.0],
            "polarization": "v",
            "yaw": 180.0
        })

    with open(args.output, "w", encoding="utf-8") as f:
        yaml.dump(elements, f, sort_keys=False)

    print(f"Generated layout with {args.rx_count} Rx and {args.tx_count} Tx elements.")
    print(f"  Centre Frequency: {args.centre_freq} GHz (lambda_0 = {lambda_0:.2f} mm)")
    print(f"  Rx Spacing: {args.rx_spacing} mm")
    print(f"  Tx Spacing: {args.tx_spacing:.2f} mm")
    print(f"Layout written to: {args.output}")
    print(f"\nTo synthesize in HFSS, run:")
    print(f"  hfss_array_builder <project_path> <source_design_name> {args.output} --centre-freq {args.centre_freq} --bandwidth {args.bandwidth}")

if __name__ == "__main__":
    main()
