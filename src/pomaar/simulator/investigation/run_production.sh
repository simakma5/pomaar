#!/usr/bin/env bash
# M4 production run of the multi-receiver bistaticity study (the investigation project, `project.PROJECT_PATH`). Run from the
# worktree root; each step logs to stdout. Order: sphere, trihedral, dihedral.
set -euo pipefail
cd "$(dirname "$0")/../../../.."
run() { ansys-vnc exec python "src/pomaar/simulator/investigation/$1.py" "${@:2}"; }
ANGLES="LIN 0.25deg 8deg 0.25deg"

# Sphere: radius 20 mm (40 mm across, like the plate targets), finely faceted (T6)
run mesh MultiRxSphere --slider 9
run sweep MultiRxSphere --sweep bistaticAngle "$ANGLES" --sweep rxAim "LIN 0 1 1" --solve
run export_sweep MultiRxSphere --fixed rxAim 0 1

# Trihedral, 40 mm edges: bistatic plane off (roll 0) and on (roll 90) its mirror plane
run sweep MultiRxTrihedral --sweep bistaticAngle "$ANGLES" --sweep rxAim "LIN 0 1 1" \
    --sweep targetRoll "LIN 0deg 90deg 90deg" --solve
run export_sweep MultiRxTrihedral --fixed rxAim 0 1 --fixed targetRoll 0deg 90deg

# Dihedral, 40 mm seam and depth: roll 0, 22.5, 45, 67.5 and 90 deg
run build_scene MultiRxDihedral --target dihedral --size 40mm
run sweep MultiRxDihedral --sweep bistaticAngle "$ANGLES" --sweep rxAim "LIN 0 1 1" \
    --sweep targetRoll "LIN 0deg 90deg 22.5deg" --solve
run export_sweep MultiRxDihedral --fixed rxAim 0 1 \
    --fixed targetRoll 0deg 22.5deg 45deg 67.5deg 90deg
echo "Production run finished."
