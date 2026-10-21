#!/usr/bin/env bash
# Cancels AEDT's 'file is in use' dialogs (401 x 119 px, Cancel at 1282, 764 on this display) until
# no Python process runs in the container. Needs the project's x11 extra (python-xlib):
#   uv sync --extra x11
cd "$(dirname "$0")/../../../../.."
X="ansys-vnc exec python src/pomaar/simulator/investigation/validation"
n=0
while ansys-vnc ps | grep -q " python "; do
    if $X/x11_windows.py 2>/dev/null | grep -q "401x119"; then
        $X/x11_click.py 1282 764 >/dev/null 2>&1
        n=$((n + 1)); echo "cancelled dialog $n"
    fi
    sleep 5
done
echo "no Python process left; $n dialogs cancelled"
