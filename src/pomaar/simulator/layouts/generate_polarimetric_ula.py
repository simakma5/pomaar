#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import yaml

def main():
    parser = argparse.ArgumentParser(description="Generate a Polarimetric ULA MIMO Layout YAML.")
    parser.add_argument("--rx-count", type=int, default=4, help="Number of Rx elements (default: 4)")
    parser.add_argument("--tx-count", type=int, default=4, help="Number of Tx elements (default: 4)")
    parser.add_argument("--center-freq", "--centre-freq", type=float, default=79.0, help="Center frequency in GHz (default: 79.0)")
    parser.add_argument("--bandwidth", type=float, default=4.0, help="Sweep bandwidth in GHz (default: 4.0)")
    parser.add_argument("--rx-spacing", type=float, default=None, help="Spacing between Rx elements in mm (default: 0.5 * lambda_0)")
    parser.add_argument("--tx-spacing", type=float, default=None, help="Spacing between Tx elements in mm (default: rx_count * rx_spacing)")
    parser.add_argument("--tx-offset-y", type=float, default=12.0, help="Y-offset of the Tx array in mm (default: 12.0)")
    parser.add_argument("-o", "--output", default="polarimetric_ula.yaml", help="Output YAML filename (default: polarimetric_ula.yaml)")

    args = parser.parse_args()

    # Dynamically calculate spacing based on center frequency if not specified
    lambda_0 = 299.792458 / args.center_freq
    if args.rx_spacing is None:
        args.rx_spacing = round(0.5 * lambda_0, 2)

    if args.tx_spacing is None:
        args.tx_spacing = args.rx_count * args.rx_spacing

    elements = []
    # Dual-pol Receivers (add V and H elements at each location)
    rx_offset_x = (args.rx_count - 1) * args.rx_spacing / 2.0
    for i in range(args.rx_count):
        x_pos = i * args.rx_spacing - rx_offset_x
        rx_idx_offset = f"({i} - ({args.rx_count} - 1) / 2.0)" if args.rx_count > 1 else "0"
        elements.append({
            "label": f"Rx_{i+1}",
            "role": "Rx",
            "position": [round(x_pos, 4), 0.0, 0.0],
            "position_expression": [f"{rx_idx_offset} * rxSpacing", "0mm", "0mm"],
            "polarization": "v",
            "yaw": 0.0
        })
        elements.append({
            "label": f"Rx_{i+1}",
            "role": "Rx",
            "position": [round(x_pos, 4), 0.0, 0.0],
            "position_expression": [f"{rx_idx_offset} * rxSpacing", "0mm", "0mm"],
            "polarization": "h",
            "yaw": 90.0
        })

    # Alternating polarizations for Transmitters
    tx_offset_x = (args.tx_count - 1) * args.tx_spacing / 2.0
    for i in range(args.tx_count):
        x_pos = i * args.tx_spacing - tx_offset_x
        tx_idx_offset = f"({i} - ({args.tx_count} - 1) / 2.0)" if args.tx_count > 1 else "0"
        is_v_pol = (i % 2 == 0)
        elements.append({
            "label": f"Tx_{i+1}",
            "role": "Tx",
            "position": [round(x_pos, 4), round(args.tx_offset_y, 4), 0.0],
            "position_expression": [f"{tx_idx_offset} * txSpacing", "txOffsetY", "0mm"],
            "polarization": "v" if is_v_pol else "h",
            "yaw": 180.0 if is_v_pol else 270.0
        })

    layout_data = {
        "metadata": {
            "topology": "polarimetric_ula",
            "center_frequency_ghz": args.center_freq,
            "bandwidth_ghz": args.bandwidth,
            "variables": {
                "rxSpacing": f"{args.rx_spacing:.4f}mm",
                "txSpacing": f"{args.tx_spacing:.4f}mm",
                "txOffsetY": f"{args.tx_offset_y:.4f}mm",
                "pcbMargin": "0.0mm",
            },
            "board": {
                "width_formula": f"max(({args.rx_count} - 1) * rxSpacing, ({args.tx_count} - 1) * txSpacing) + 2 * (unitCellExtentX + pcbMargin)",
                "length_formula": "txOffsetY + 2 * (unitCellExtentY + pcbMargin)",
            },
        },
        "elements": elements,
    }

    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = Path(__file__).resolve().parent / output_path

    with open(output_path, "w", encoding="utf-8") as f:
        yaml.dump(layout_data, f, sort_keys=False)

    print(f"Generated layout with {args.rx_count} Rx (dual-pol) and {args.tx_count} Tx (alternating V/H) elements.")
    print(f"  Center Frequency: {args.center_freq} GHz (lambda_0 = {lambda_0:.2f} mm)")
    print(f"  Rx Spacing (Dual-Pol): {args.rx_spacing} mm")
    print(f"  Tx Spacing (Alternating V/H): {args.tx_spacing:.2f} mm")
    print(f"Layout written to: {output_path}")
    print(f"\nTo synthesize in HFSS, run:")
    print(f"  hfss_array_builder <project_path> <source_design_name> {output_path} --center-freq {args.center_freq} --bandwidth {args.bandwidth}")

if __name__ == "__main__":
    main()
