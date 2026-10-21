"""Saves a PNG of the region (x, y, w, h) of display :1.

    ansys-vnc exec python src/pomaar/agent/x11_shot.py X Y W H out.png

Needs the project's x11 extra (`uv sync --extra x11`).
"""

import sys

from PIL import Image
from Xlib import X, display

x, y, w, h = (int(v) for v in sys.argv[1:5])
out = sys.argv[5]
root = display.Display(":1").screen().root
raw = root.get_image(x, y, w, h, X.ZPixmap, 0xFFFFFFFF)
Image.frombytes("RGB", (w, h), raw.data, "raw", "BGRX").save(out)
print("saved", out)
