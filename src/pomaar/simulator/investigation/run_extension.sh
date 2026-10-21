#!/usr/bin/env bash
# Extension runs of the bistaticity study (sandbox project only), after the supervisor's
# suggestions: more reflectors (slant trihedral, flat plate) and the differential reflectivity
# from split antenna pairs. Receivers aimed (rxAim 1) only: aiming changes the span only
# (reference/rx_aim). Run from the repository root; each step logs to stdout.
set -euo pipefail
cd "$(dirname "$0")/../../../.."
run() { ansys-vnc exec python "src/pomaar/simulator/investigation/$1.py" "${@:2}"; }
ANGLES="LIN 0.25deg 8deg 0.25deg"
PAIR_ANGLES="LIN 0.5deg 6deg 0.5deg"

# 1. Slant trihedral: roll 45 deg, between its vertical (0) and horizontal (90) mirror cases
run sweep MultiRxTrihedral --name Roll45Sweep --sweep bistaticAngle "$ANGLES" \
    --sweep rxAim 1 --sweep targetRoll 45deg --solve
run export_sweep MultiRxTrihedral --fixed rxAim 1 --fixed targetRoll 45deg

# 2. Flat square plate, 40 mm edges, facing the radar: single bounce with a lobe
run build_scene MultiRxPlate --target plate --size 40mm
run sweep MultiRxPlate --sweep bistaticAngle "$ANGLES" --solve
run export_sweep MultiRxPlate --fixed rxAim 1

# 3. Split pairs (TxV-RxV, gap, TxH-RxH) and the co-located reference pair; gap 0 and 2 b
#    (pair bisectors beta and 3 beta apart)
run build_scene PairsSphere --layout pairs --target sphere --size 20mm
run mesh PairsSphere --slider 9
run sweep PairsSphere --sweep pairAngle "$PAIR_ANGLES" --sweep gapRatio "LIN 0 2 2" --solve
run export_sweep PairsSphere --sweep-variable pairAngle --fixed gapRatio 0 2

run build_scene PairsTrihedral --layout pairs --target trihedral --size 40mm
run sweep PairsTrihedral --sweep pairAngle "$PAIR_ANGLES" --sweep gapRatio "LIN 0 2 2" \
    --sweep targetRoll "LIN 0deg 90deg 90deg" --solve
run export_sweep PairsTrihedral --sweep-variable pairAngle --fixed gapRatio 0 2 \
    --fixed targetRoll 0deg 90deg

run build_scene PairsDihedral --layout pairs --target dihedral --size 40mm
run sweep PairsDihedral --sweep pairAngle "$PAIR_ANGLES" --sweep gapRatio "LIN 0 2 2" \
    --sweep targetRoll "LIN 0deg 22.5deg 22.5deg" --solve
run export_sweep PairsDihedral --sweep-variable pairAngle --fixed gapRatio 0 2 \
    --fixed targetRoll 0deg 22.5deg
echo "Extension run finished."
