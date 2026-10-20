#!/usr/bin/env python3
"""
Exports the solved `bistaticAngle` sweep of a sandbox design once per combination of the other
swept variables (rxAim, targetRoll, ...), which are set as nominal values first:

    ansys-vnc exec python src/pomaar/simulator/investigation/export_sweep.py MultiRxDihedral \
        --fixed rxAim 0 1 --fixed targetRoll 0deg 45deg

Each export (<sandbox>/exports/<design>_<var><value>...npz) carries the scene geometry in its
metadata, in the ISO frame of the scene (x boresight, y left, z up), per variation:
`positions_m` {port: [[x, y, z], ...]} and `target_m` [[x, y, z], ...] (see `scene_geometry`).
"""

import argparse
import itertools
import json
import re

import numpy as np

from pomaar.simulator.array_templates import parse_port
from pomaar.simulator.investigation.sandbox import SANDBOX_DIR, open_sandbox, release
from pomaar.simulator.sbr_results import export_sbr_parametric_sweep

UNITS = {"meter": 1.0, "m": 1.0, "mm": 1e-3, "cm": 1e-2, "deg": 1.0, "": 1.0}


def value(text):
    """'4meter' -> 4.0, '20mm' -> 0.02 (lengths in m), '1deg' -> 1.0 (angles in deg)."""
    match = re.fullmatch(r"\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*([a-z]*)\s*", str(text))
    if not match or match[2] not in UNITS:
        raise ValueError(f"Cannot read the AEDT value '{text}'.")
    return float(match[1]) * UNITS[match[2]]


def scene_geometry(port, variables):
    """
    Phase-centre position (m, ISO frame) of a port of `build_scene`'s scene: the Tx site at the
    origin, Rx site k at y = targetRange tan(k bistaticAngle), H at z = +spacingHV/2 and V at
    -spacingHV/2.
    """
    role, pol, key = parse_port(port)
    k = int(key.strip("_").split("_")[1]) if role == "Rx" else 0
    y = value(variables["targetRange"]) * np.tan(np.radians(k * value(variables["bistaticAngle"])))
    z = (0.5 if pol == "H" else -0.5) * value(variables["spacingHV"])
    return [0.0, y, z]


def target_position(variables):
    """TargetCS origin (m, ISO frame): range, azimuth from +x towards +y, elevation."""
    r = value(variables["targetRange"])
    az, el = (
        np.radians(value(variables["targetAzimuth"])),
        np.radians(value(variables["targetElevation"])),
    )
    return [r * np.cos(el) * np.cos(az), r * np.cos(el) * np.sin(az), r * np.sin(el)]


def attach_geometry(path, ports):
    data = dict(np.load(path))
    metadata = json.loads(str(data["metadata"]))
    nominal = dict(metadata["nominal_variables"])
    positions = {port: [] for port in ports}
    targets = []
    for angle in data["sweep_values"]:
        variables = {**nominal, metadata["sweep_variable"]: f"{angle}deg"}
        for port in ports:
            positions[port].append(scene_geometry(port, variables))
        targets.append(target_position(variables))
    metadata.update(frame="ISO 8855 (x boresight)", positions_m=positions, target_m=targets)
    data["metadata"] = json.dumps(metadata)
    np.savez(path, **data)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("design")
    parser.add_argument("--fixed", nargs="+", action="append", default=[], metavar="VAR VALUES")
    parser.add_argument("--sweep-variable", default="bistaticAngle")
    args = parser.parse_args(argv)

    names = [f[0] for f in args.fixed]
    combinations = list(itertools.product(*[f[1:] for f in args.fixed])) or [()]
    desktop, app = open_sandbox(args.design)
    try:
        for combination in combinations:
            for name, setting in zip(names, combination):
                app[name] = setting
            app.save_project()
            label = "_".join(
                f"{name}{re.sub(r'[^0-9.]', '', setting)}"
                for name, setting in zip(names, combination)
            )
            path = SANDBOX_DIR / "exports" / f"{args.design}{'_' + label if label else ''}.npz"
            _, _, _, expressions = export_sbr_parametric_sweep(
                app, str(path), sweep_variable=args.sweep_variable
            )
            ports = sorted({p for e in expressions for p in re.findall(r"[\w:]+", e)[1:]})
            attach_geometry(path, ports)
            print(f"  + geometry of {len(ports)} ports")
    finally:
        release(desktop)


if __name__ == "__main__":
    main()
