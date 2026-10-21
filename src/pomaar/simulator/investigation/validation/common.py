"""
Guard of the pipeline validation: scripts attach to the live AEDT GUI session and open only the
validation copy of the project
(~/Projects/AEDT/mimo-polarimetry-validation/mimo_polarimetry_validation.aedt), never the
original. The harness's logs and exports go to data/interim/validation/ (`VALIDATION_DIR`).
"""

from pathlib import Path

from ansys.aedt.core import Hfss

from pomaar.config import INTERIM_DATA_DIR
from pomaar.simulator.hfss_array_builder import connect_live_session

VALIDATION_DIR = INTERIM_DATA_DIR / "validation"
PROJECT_NAME = "mimo_polarimetry_validation"
PROJECT_PATH = Path.home() / "Projects/AEDT/mimo-polarimetry-validation" / f"{PROJECT_NAME}.aedt"
ORIGINALS = [Path.home() / "Projects/AEDT/mimo-polarimetry/mimo_polarimetry.aedt"]


def open_validation(design=None):
    """(desktop, app) for `design` of the validation project, opened in the GUI if needed."""
    if any(PROJECT_PATH.resolve() == o.resolve() for o in ORIGINALS):
        raise RuntimeError("The validation path points at an original project.")
    desktop = connect_live_session()
    if PROJECT_NAME not in desktop.project_list:
        desktop.odesktop.OpenProject(str(PROJECT_PATH))
    app = Hfss(project=PROJECT_NAME, design=design, new_desktop=False, close_on_exit=False)
    path = Path(app.project_path) / f"{app.project_name}.aedt"
    if path.resolve() != PROJECT_PATH.resolve():
        raise RuntimeError(f"Refusing to work on {path}.")
    return desktop, app


def release(desktop):
    desktop.release_desktop(close_projects=False, close_on_exit=False)
