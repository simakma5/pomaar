"""Lists viewable windows on :1 between 200x60 and 900x400 px with position and WM_CLASS."""

from Xlib import X, display

d = display.Display(":1")
root = d.screen().root


def windows(w, depth=0):
    for child in w.query_tree().children:
        yield child, depth
        yield from windows(child, depth + 1)


for w, depth in windows(root):
    try:
        attrs = w.get_attributes()
        geom = w.get_geometry()
        pos = w.translate_coords(root, 0, 0)
        cls = w.get_wm_class()
    except Exception:
        continue
    if attrs.map_state == X.IsViewable and 200 <= geom.width <= 900 and 60 <= geom.height <= 400:
        print(f"{w.id:#x} depth {depth} {geom.width}x{geom.height} at ({-pos.x},{-pos.y}) {cls}")
