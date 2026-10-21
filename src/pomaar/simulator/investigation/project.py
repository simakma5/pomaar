"""
The AEDT project of the bistaticity investigation and the folders of its products.

* `PROJECT_PATH`: the parametric-beam scenes,
  ~/Projects/AEDT/bistaticity-parametric-beams/bistaticity_parametric_beams.aedt (a copy of
  ~/Projects/AEDT/bistaticity-analysis/bistaticity_analysis.aedt with the investigation's
  designs). The AEDT-side scripts attach to the live GUI session through `open_project`, which
  opens the project there if needed; their `--project` (`add_project_argument`) takes another.
* <exports>, `EXPORTS_DIR`: data/interim/bistaticity/parametric_beams, the .npz exports and,
  in logs/, the run logs.
* <results>, `RESULTS_DIR`: data/processed/bistaticity/parametric_beams, the summaries, the
  split pairs, the scaling and the receiver-aiming check; the per-export analyses
  (`bistaticity_analysis --per-receiver`) go to its per_export/ (`ANALYSIS_DIR`).

The scenes use idealized parametric beams (no embedded patterns, no antenna cross-pol), hence
the parametric_beams layer next to the full-wave cluster results in data/*/bistaticity/.
"""

from pathlib import Path

from ansys.aedt.core import Hfss

from pomaar.config import INTERIM_DATA_DIR, PROCESSED_DATA_DIR
from pomaar.simulator.hfss_array_builder import connect_live_session, require_designs

PROJECT_PATH = (
    Path.home() / "Projects/AEDT/bistaticity-parametric-beams/bistaticity_parametric_beams.aedt"
)
EXPORTS_DIR = INTERIM_DATA_DIR / "bistaticity" / "parametric_beams"
RESULTS_DIR = PROCESSED_DATA_DIR / "bistaticity" / "parametric_beams"
ANALYSIS_DIR = RESULTS_DIR / "per_export"


def add_project_argument(parser):
    parser.add_argument(
        "--project",
        type=Path,
        default=PROJECT_PATH,
        help=f"AEDT project (default: {PROJECT_PATH})",
    )


def open_project(design=None, project=PROJECT_PATH):
    """
    Attaches to the GUI session, opens `project` there if needed and returns (desktop, app) for
    `design` (an existing design; None for the active one).
    """
    project = Path(project).expanduser().resolve()
    if not project.is_file():
        raise FileNotFoundError(f"No AEDT project at {project}.")
    desktop = connect_live_session()
    if project.stem not in desktop.project_list:
        print(f"Opening {project} in the GUI session...")
        desktop.odesktop.OpenProject(str(project))
    if design is not None:
        require_designs(desktop, project.stem, [design])
    app = Hfss(project=project.stem, design=design, new_desktop=False, close_on_exit=False)
    opened = Path(app.project_path) / f"{app.project_name}.aedt"
    if opened.resolve() != project:
        raise RuntimeError(f"'{project.stem}' is open from {opened}, not from {project}.")
    return desktop, app


def reopen_project(desktop, app, design=None):
    """
    Saves, closes and reopens the project of `app`, which releases AEDT's handles on the result
    files of deleted designs (a design re-created under a deleted design's name otherwise meets
    a modal 'file ... is in use' dialog on save); returns a fresh app for `design`.
    """
    name, path = app.project_name, Path(app.project_path) / f"{app.project_name}.aedt"
    desktop.odesktop.SetActiveProject(name).Save()
    desktop.odesktop.CloseProject(name)
    desktop.odesktop.OpenProject(str(path))
    return Hfss(project=name, design=design, new_desktop=False, close_on_exit=False)


def release(desktop):
    """Releases the session connection; AEDT and its projects stay open."""
    desktop.release_desktop(close_projects=False, close_on_exit=False)
