#!/usr/bin/env python3
"""
Generates a 4-element polarimetric quad-cluster layout (Tx_H, Tx_V, Rx_H, Rx_V)
with parameterized radial distance from the geometrical centre.
"""

import argparse
import math
from pathlib import Path
import sys
import yaml


def main():
    parser = argparse.ArgumentParser(
        description="Generate a 4-Element Polarimetric MIMO Quad-Cluster Layout YAML (Tx_H, Tx_V, Rx_H, Rx_V)."
    )
    parser.add_argument(
        "--radius",
        "-r",
        type=float,
        default=None,
        help="Radial distance of elements from the center in mm (default: 0.35 * lambda_0 => copol baseline = 0.7 * lambda_0)",
    )
    parser.add_argument(
        "--radius-lambda",
        type=float,
        default=None,
        help="Radial distance in units of wavelengths (lambda_0) (default: 0.35)",
    )
    parser.add_argument(
        "--copol-spacing",
        type=float,
        default=None,
        help="Co-polar baseline spacing (2*r) in mm (default: 0.7 * lambda_0)",
    )
    parser.add_argument(
        "--center-freq",
        "--centre-freq",
        type=float,
        default=79.0,
        help="Center frequency in GHz (default: 79.0)",
    )
    parser.add_argument(
        "--bandwidth",
        type=float,
        default=4.0,
        help="Sweep bandwidth in GHz (default: 4.0)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="polarimetric_cluster.yaml",
        help="Output YAML filename (default: polarimetric_cluster.yaml)",
    )

    args = parser.parse_args()

    # Speed of light in vacuum (mm/s) -> wavelength in mm
    lambda_0 = 299.792458 / args.center_freq

    # Determine radial distance (mm)
    if args.radius is not None:
        radius_mm = args.radius
    elif args.copol_spacing is not None:
        radius_mm = args.copol_spacing / 2.0
    elif args.radius_lambda is not None:
        radius_mm = args.radius_lambda * lambda_0
    else:
        # Default: co-polar baseline = 0.7 * lambda_0 => radius = 0.35 * lambda_0
        radius_mm = round(0.35 * lambda_0, 2)

    copol_spacing_mm = 2.0 * radius_mm
    crosspol_spacing_mm = math.sqrt(2.0) * radius_mm

    # Build 4-element cluster with opposed co-polar pairs:
    # - Horizontal pair (H-pol) along X-axis:
    #   * Tx_H at (+r, 0): yaw = 90 deg (feedline to +X, H-pol)
    #   * Rx_H at (-r, 0): yaw = 270 deg (feedline to -X, H-pol)
    #   * Co-polar HH baseline: d_co = 2*r along X, virtual phase center = [0, 0, 0]
    # - Vertical pair (V-pol) along Y-axis:
    #   * Tx_V at (0, +r): yaw = 180 deg (feedline to +Y, V-pol)
    #   * Rx_V at (0, -r): yaw = 0 deg (feedline to -Y, V-pol)
    #   * Co-polar VV baseline: d_co = 2*r along Y, virtual phase center = [0, 0, 0]
    layout_data = {
        "metadata": {
            "topology": "polarimetric_cluster",
            "center_frequency_ghz": args.center_freq,
            "bandwidth_ghz": args.bandwidth,
            "variables": {
                "copolarSpacing": f"{copol_spacing_mm:.4f}mm",
                "clusterRadius": "copolarSpacing / 2",
                "crosspolarSpacing": "sqrt(2) * clusterRadius",
                "pcbMargin": "0.0mm",
            },
            "board": {
                "width_formula": "2 * (clusterRadius + unitCellExtent + pcbMargin)",
                "length_formula": "2 * (clusterRadius + unitCellExtent + pcbMargin)",
            },
        },
        "elements": [
            {
                "label": "Tx",
                "role": "Tx",
                "polarization": "h",
                "position": [round(radius_mm, 4), 0.0, 0.0],
                "position_expression": ["clusterRadius", "0mm", "0mm"],
                "yaw": 90.0,
            },
            {
                "label": "Rx",
                "role": "Rx",
                "polarization": "h",
                "position": [round(-radius_mm, 4), 0.0, 0.0],
                "position_expression": ["-clusterRadius", "0mm", "0mm"],
                "yaw": 270.0,
            },
            {
                "label": "Tx",
                "role": "Tx",
                "polarization": "v",
                "position": [0.0, round(radius_mm, 4), 0.0],
                "position_expression": ["0mm", "clusterRadius", "0mm"],
                "yaw": 180.0,
            },
            {
                "label": "Rx",
                "role": "Rx",
                "polarization": "v",
                "position": [0.0, round(-radius_mm, 4), 0.0],
                "position_expression": ["0mm", "-clusterRadius", "0mm"],
                "yaw": 0.0,
            },
        ],
    }

    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = Path(__file__).resolve().parent / output_path

    with open(output_path, "w", encoding="utf-8") as f:
        yaml.dump(layout_data, f, sort_keys=False)

    print(
        f"Generated 4-element polarimetric cluster layout at {args.center_freq} GHz (lambda_0 = {lambda_0:.2f} mm):"
    )
    print(f"  Radial Distance (r):        {radius_mm:.3f} mm ({radius_mm / lambda_0:.3f} * lambda_0)")
    print(
        f"  Co-Polar Baseline (d_co):   {copol_spacing_mm:.3f} mm ({copol_spacing_mm / lambda_0:.3f} * lambda_0)"
    )
    print(
        f"  Cross-Polar Baseline (d_x): {crosspol_spacing_mm:.3f} mm ({crosspol_spacing_mm / lambda_0:.3f} * lambda_0)"
    )
    print(f"  Virtual Phase Centers:      r_v(HH) = r_v(VV) = [0.0, 0.0, 0.0]")
    print(f"Layout written to: {output_path}")
    print("\nTo synthesize in HFSS, run:")
    print(
        f"  hfss_array_builder <project_path> <source_design_name> {output_path} --center-freq {args.center_freq} --bandwidth {args.bandwidth}"
    )


if __name__ == "__main__":
    main()
