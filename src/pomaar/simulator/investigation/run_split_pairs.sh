#!/usr/bin/env bash
# Runs of 2026-10-28/29 (the investigation project, `project.PROJECT_PATH`): split-pair scenes at
# gap 0 for the dihedral and trihedral at roll 45 deg and the plate (`pairs.py`), and
# receiver-on-boresight ladders (rxAim 0) of the plate and the trihedral at roll 45, 270 and
# 315 deg. The 270/315 deg ladders served only a split-pair reconstruction that was dropped
# (investigation_log.md, 2026-10-29). New parametric setups only; the existing ones and their
# results stay. Run from the repository root; each step logs to stdout.
set -euo pipefail
cd "$(dirname "$0")/../../../.."
run() { ansys-vnc exec python "src/pomaar/simulator/investigation/$1.py" "${@:2}"; }
ANGLES="LIN 0.25deg 8deg 0.25deg"
PAIR_ANGLES="LIN 0.5deg 6deg 0.5deg"

# 1. Split pairs at gap 0: dihedral and trihedral at roll 45 deg, the plate
run sweep PairsDihedral --name Roll45Sweep --sweep pairAngle "$PAIR_ANGLES" --sweep gapRatio 0 \
    --sweep targetRoll 45deg --solve
run export_sweep PairsDihedral --sweep-variable pairAngle --fixed gapRatio 0 --fixed targetRoll 45deg

run sweep PairsTrihedral --name Roll45Sweep --sweep pairAngle "$PAIR_ANGLES" --sweep gapRatio 0 \
    --sweep targetRoll 45deg --solve
run export_sweep PairsTrihedral --sweep-variable pairAngle --fixed gapRatio 0 --fixed targetRoll 45deg

run build_scene PairsPlate --layout pairs --target plate --size 40mm
run sweep PairsPlate --sweep pairAngle "$PAIR_ANGLES" --sweep gapRatio 0 --solve
run export_sweep PairsPlate --sweep-variable pairAngle --fixed gapRatio 0

# 2. Receivers on boresight (rxAim 0): the plate, the trihedral at 45 deg and at 270 and 315 deg
run sweep MultiRxPlate --name Rx0Sweep --sweep bistaticAngle "$ANGLES" --sweep rxAim 0 --solve
run export_sweep MultiRxPlate --fixed rxAim 0
for roll in 45 270 315; do
    run sweep MultiRxTrihedral --name "Rx0Roll${roll}Sweep" --sweep bistaticAngle "$ANGLES" \
        --sweep rxAim 0 --sweep targetRoll "${roll}deg" --solve
done
run export_sweep MultiRxTrihedral --fixed rxAim 0 --fixed targetRoll 45deg 270deg 315deg
echo "Split-pair run finished."
