# Agent tools for the live AEDT session

Small tools for driving the AEDT 2025 R2 GUI session in the `ansys-vnc` container, written
for a coding agent but just as usable by hand. Like every PyAEDT script here, they run inside
the container (`ansys-vnc exec`) and attach to the running GUI. They never start an AEDT.

| Tool | What it does | Changes anything? |
| --- | --- | --- |
| `list_projects.py` | open projects, their folders and designs | no |
| `inspect_design.py PROJECT DESIGN ...` | variables, CSs, objects, boundaries, excitations, saved setups and sweeps, reports, solutions (SBR+: scattered solution, Tx/Rx pairs); refuses missing designs instead of creating them | no |
| `x11_windows.py` | viewable windows of dialog size (200x60 to 900x400 px) on the VNC display `:1`, with position and class | no |
| `x11_shot.py X Y W H OUT.png` | screenshot of a screen region | no |
| `x11_click.py X Y` | left click at a screen position (XTEST) | yes: clicks the GUI |
| `watch_dialogs.sh` | while Python runs in the container, cancels AEDT's modal "file ... is in use" dialog (401 x 119 px, Cancel at 1282, 764 on this display) | yes: clicks Cancel |

```bash
ansys-vnc exec python -m pomaar.agent.list_projects
ansys-vnc exec python -m pomaar.agent.inspect_design ~/Projects/AEDT/<folder>/<project>.aedt Design1
ansys-vnc exec python src/pomaar/agent/x11_shot.py 0 0 1920 1080 sandbox/screen.png
```

* **Display helpers.** The X11 helpers need the project's x11 extra (`uv sync --extra x11`,
  python-xlib; Pillow comes with matplotlib).
* **The dialog.** It appears when a design is re-created under the name of one deleted in the
  same session (see the CLAUDE.md pitfalls). The scripts now avoid it by reopening the project,
  so `watch_dialogs.sh` is a fallback. Its coordinates are for this display: check them with
  `x11_windows.py` and `x11_shot.py` first.
* **Scratch output.** Screenshots, logs and throwaway exports go to `sandbox/` at the repository
  root (gitignored).
* **Project copies.** AEDT project copies for experiments go to `~/Projects/AEDT/` under a
  telling name, never into the repository.
