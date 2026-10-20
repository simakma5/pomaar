#!/usr/bin/env bash
# Size and range scaling of the multi-receiver results (sandbox project only): trihedral at roll
# 0 (first-order cross-pol) at 80 mm and 160 mm (4 m), 80 mm at 16 m (far field, 2 D^2 / lambda
# = 6.7 m), 40 mm at 1 m (near field), and a 160 mm dihedral at roll 90 (retro plane). rxAim 1.
set -euo pipefail
cd "$(dirname "$0")/../../../.."
run() { ansys-vnc exec python "src/pomaar/simulator/investigation/$1.py" "${@:2}"; }

scene() {  # design, target, size, range, roll, angles ('-' as target: design exists)
    [ "$2" = - ] || run build_scene "$1" --target "$2" --size "$3" --range "$4" --roll "$5"
    run sweep "$1" --sweep bistaticAngle "$6" --solve
    run export_sweep "$1" --fixed rxAim 1
}
# MultiRxTrihedralS160 exists already (built for the timing solve)
scene MultiRxTrihedralS160 - 160mm 4m 0deg "LIN 0.25deg 3deg 0.25deg"
scene MultiRxTrihedralS80 trihedral 80mm 4m 0deg "LIN 0.25deg 3deg 0.25deg"
scene MultiRxTrihedralS80R16 trihedral 80mm 16m 0deg "LIN 0.25deg 3deg 0.25deg"
scene MultiRxTrihedralR1 trihedral 40mm 1m 0deg "LIN 0.5deg 6deg 0.5deg"
scene MultiRxDihedralS160 dihedral 160mm 4m 90deg "LIN 0.25deg 3deg 0.25deg"
echo "Scaling run finished."
