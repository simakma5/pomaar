"""Clicks the left mouse button at screen position (x, y) on :1 through XTEST.

    ansys-vnc exec python src/pomaar/agent/x11_click.py X Y

Needs the project's x11 extra (`uv sync --extra x11`).
"""

import sys
import time

from Xlib import X, display
from Xlib.ext import xtest

d = display.Display(":1")
x, y = int(sys.argv[1]), int(sys.argv[2])
xtest.fake_input(d, X.MotionNotify, x=x, y=y)
d.sync()
time.sleep(0.2)
xtest.fake_input(d, X.ButtonPress, 1)
xtest.fake_input(d, X.ButtonRelease, 1)
d.sync()
print("clicked", x, y)
