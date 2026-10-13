"""
Unit-cell binding and output naming shared by the layout generators.

A layout is bound to the unit-cell designs it is built from: `metadata.primary_unit_cell` names
the design every element uses unless `metadata.unit_cells` ({h, v, hv: design}) or an element's
own `unit_cell` names another, and `hfss_array_builder <project.aedt> <layout.yaml>` takes them
from there. The file is named after its primary unit cell, `<UnitCell>_<topology>.yaml`, so one
generator run per unit cell gives one layout per cell, e.g.

    for cell in ApertureCoupledPatch SeriesFedPatch; do
        python generate_polarimetric_cluster.py --unit-cell $cell
        ansys-vnc exec hfss_array_builder <project.aedt> ${cell}_polarimetric_cluster.yaml
    done

and the builder names the array design <UnitCell><Topology>, e.g.
ApertureCoupledPatchPolarimetricCluster.
The frequency sweep is not part of a layout: the builder copies the primary unit cell's setup.
"""

from pathlib import Path

import yaml


def add_unit_cell_arguments(parser):
    """
    Adds --unit-cell (the primary design, required), --extra-unit-cell POL=DESIGN,
    --drop-unused-ports, --phase-reference-mode and -o/--output.
    """
    parser.add_argument(
        "--unit-cell",
        required=True,
        metavar="DESIGN",
        help=(
            "Primary unit-cell design of the layout (metadata.primary_unit_cell); it also names "
            "the array design <DESIGN><LayoutFileStem>"
        ),
    )
    parser.add_argument(
        "--extra-unit-cell",
        action="append",
        default=[],
        metavar="POL=DESIGN",
        help=(
            "Further unit cell for h, v or hv (dual-polarized) elements, e.g. hv=DualPatch "
            "(metadata.unit_cells; repeatable)"
        ),
    )
    parser.add_argument(
        "--drop-unused-ports",
        action="store_true",
        default=False,
        help=(
            "Remove the unused ports of dual-port unit cells instead of terminating them in "
            "matched '_Load' wave ports (metadata.drop_unused_ports)"
        ),
    )
    parser.add_argument(
        "--phase-reference-mode",
        choices=["per_port", "single"],
        default="per_port",
        help=(
            "Far-field phase reference for SBR+ linking: each port at its own phase centre "
            "(per_port, default) or all at the array centroid (single; far-field targets only)"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output YAML file (default: <UnitCell>_<topology>.yaml next to the generator)",
    )


def unit_cell_metadata(parser, args):
    """The layout metadata entries for the unit cells chosen on the command line."""
    metadata = {"primary_unit_cell": args.unit_cell}
    extra = {}
    for value in args.extra_unit_cell:
        key, sep, design = value.partition("=")
        key = key.strip().lower()
        if not sep or not design.strip() or key not in ("h", "v", "hv", "vh"):
            parser.error(f"--extra-unit-cell expects h|v|hv=<Design>, got '{value}'")
        extra["hv" if key == "vh" else key] = design.strip()
    if extra:
        metadata["unit_cells"] = extra
    if args.drop_unused_ports:
        metadata["drop_unused_ports"] = True
    if args.phase_reference_mode != "per_port":
        metadata["phase_reference_mode"] = args.phase_reference_mode
    return metadata


def layout_filename(args, topology, tag=""):
    """The output file: -o as given, else '<UnitCell>_<topology><tag>.yaml'."""
    if args.output:
        path = Path(args.output)
        if tag and not path.stem.endswith(tag):
            path = path.with_name(f"{path.stem}{tag}{path.suffix}")
        return path
    return Path(f"{args.unit_cell}_{topology}{tag}.yaml")


def write_layout(layout_data, output, script_file):
    """Writes the layout YAML (relative paths next to the generator) and returns its path."""
    output_path = Path(output)
    if not output_path.is_absolute():
        output_path = Path(script_file).resolve().parent / output_path
    with open(output_path, "w", encoding="utf-8") as f:
        yaml.dump(layout_data, f, sort_keys=False)
    print(f"Layout written to: {output_path}")
    print("\nTo synthesize in HFSS (inside the AEDT container, attached to the GUI session), run:")
    print(f"  hfss_array_builder <project.aedt> {output_path}")
    return output_path
