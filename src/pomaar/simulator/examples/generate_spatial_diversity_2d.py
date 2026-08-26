#!/usr/bin/env python3
import argparse
import sys
import yaml

def main():
    parser = argparse.ArgumentParser(description="Generate a 2D Spatial Diversity MIMO Layout YAML.")
    parser.add_argument("--rx-count", type=int, default=3, help="Number of Rx elements per dimension in the L-shape (default: 3)")
    parser.add_argument("--tx-count", type=int, default=3, help="Number of Tx elements (default: 3)")
    parser.add_argument("--center-freq", "--centre-freq", type=float, default=79.0, help="Center frequency in GHz (default: 79.0)")
    parser.add_argument("--bandwidth", type=float, default=4.0, help="Sweep bandwidth in GHz (default: 4.0)")
    parser.add_argument("--rx-spacing-x", type=float, default=None, help="Rx spacing along X in mm (default: 0.5 * lambda_0)")
    parser.add_argument("--rx-spacing-y", type=float, default=None, help="Rx spacing along Y in mm (default: 0.5 * lambda_0)")
    parser.add_argument("--tx-spacing-x", type=float, default=None, help="Tx spacing along X in mm (default: 1.5 * lambda_0)")
    parser.add_argument("--tx-spacing-y", type=float, default=6.0, help="Tx spacing offset in Y in mm (default: 6.0)")
    parser.add_argument("--tx-offset-y", type=float, default=12.0, help="Base Y-offset of the Tx array in mm (default: 12.0)")
    parser.add_argument("-o", "--output", default="spatial_diversity_2d_layout.yaml", help="Output YAML filename (default: spatial_diversity_2d_layout.yaml)")

    args = parser.parse_args()

    # Dynamically calculate spacing based on center frequency if not specified
    lambda_0 = 299.792458 / args.center_freq
    if args.rx_spacing_x is None:
        args.rx_spacing_x = round(0.5 * lambda_0, 2)
    if args.rx_spacing_y is None:
        args.rx_spacing_y = round(0.5 * lambda_0, 2)
    if args.tx_spacing_x is None:
        args.tx_spacing_x = round(1.5 * lambda_0, 2)

    # Dynamic L-shaped Rx Array
    elements = []
    
    # Azimuth elements along X
    rx_offset_x = (args.rx_count - 1) * args.rx_spacing_x / 2.0
    for i in range(args.rx_count):
        x_pos = i * args.rx_spacing_x - rx_offset_x
        elements.append({
            "label": f"Rx_Az_{i+1}",
            "role": "Rx",
            "pos": [x_pos, 0.0, 0.0],
            "polarization": "v",
            "yaw": 0.0
        })

    # Elevation elements along Y (starting from y = rx_spacing_y, positioned at x = 0.0)
    for j in range(1, args.rx_count):
        elements.append({
            "label": f"Rx_El_{j}",
            "role": "Rx",
            "pos": [0.0, j * args.rx_spacing_y, 0.0],
            "polarization": "v",
            "yaw": 0.0
        })

    # Dynamic Tx Array
    tx_offset_x = (args.tx_count - 1) * args.tx_spacing_x / 2.0 if args.tx_count > 1 else 0.0
    for i in range(args.tx_count):
        x_pos = i * args.tx_spacing_x - tx_offset_x if args.tx_count > 1 else 0.0
        y_pos = args.tx_offset_y
        if args.tx_count > 2 and i == args.tx_count // 2:
            y_pos += args.tx_spacing_y
            
        elements.append({
            "label": f"Tx_{i+1}",
            "role": "Tx",
            "pos": [x_pos, y_pos, 0.0],
            "polarization": "v",
            "yaw": 180.0
        })

    with open(args.output, "w", encoding="utf-8") as f:
        yaml.dump(elements, f, sort_keys=False)

    print(f"Generated 2D layout with {args.rx_count} Rx (L-shape) and {args.tx_count} Tx (triangular/sparse) elements.")
    print(f"  Center Frequency: {args.center_freq} GHz (lambda_0 = {lambda_0:.2f} mm)")
    print(f"  Rx Spacing X: {args.rx_spacing_x} mm, Y: {args.rx_spacing_y} mm")
    print(f"  Tx Spacing X: {args.tx_spacing_x} mm, Y-Offset: {args.tx_spacing_y} mm")
    print(f"Layout written to: {args.output}")
    print(f"\nTo synthesize in HFSS, run:")
    print(f"  hfss_array_builder <project_path> <source_design_name> {args.output} --center-freq {args.center_freq} --bandwidth {args.bandwidth}")

if __name__ == "__main__":
    main()
