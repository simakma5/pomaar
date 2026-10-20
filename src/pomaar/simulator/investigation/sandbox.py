"""
Sandbox guard for the bistaticity investigation: every script attaches to the live AEDT GUI
session through `open_sandbox`, which opens only the sandbox copy of the project and refuses
anything else (in particular the original in ~/Projects/AEDT).
"""

from pathlib import Path

from ansys.aedt.core import Hfss

from pomaar.simulator.hfss_array_builder import connect_live_session

SANDBOX_DIR = Path.home() / "Repositories/pomaar/data/interim/sandbox/bistaticity"
PROJECT_NAME = "bistaticity_sandbox"
PROJECT_PATH = SANDBOX_DIR / f"{PROJECT_NAME}.aedt"
ORIGINAL_PATH = Path.home() / "Projects/AEDT/bistaticity-analysis/bistaticity_analysis.aedt"


class NotSandboxError(RuntimeError):
    """Raised when a script would act on a project other than the sandbox copy."""


def check_sandbox(app):
    """Raises NotSandboxError unless `app` is a design of the sandbox project file."""
    path = Path(app.project_path) / f"{app.project_name}.aedt"
    if app.project_name != PROJECT_NAME or path.resolve() != PROJECT_PATH.resolve():
        raise NotSandboxError(f"Refusing to work on {path}: only {PROJECT_PATH} is allowed.")


def open_sandbox(design=None):
    """
    Attaches to the GUI session, opens the sandbox project there if needed and returns
    (desktop, app) for `design` (an existing design; None for the active one).
    """
    if PROJECT_PATH.resolve() == ORIGINAL_PATH.resolve():
        raise NotSandboxError("The sandbox path points at the original project.")
    if not PROJECT_PATH.is_file():
        raise FileNotFoundError(f"No sandbox project at {PROJECT_PATH}.")
    desktop = connect_live_session()
    if PROJECT_NAME not in desktop.project_list:
        print(f"Opening {PROJECT_PATH} in the GUI session...")
        desktop.odesktop.OpenProject(str(PROJECT_PATH))
    app = Hfss(project=PROJECT_NAME, design=design, new_desktop=False, close_on_exit=False)
    check_sandbox(app)
    return desktop, app


def reopen_sandbox(desktop, design=None):
    """
    Saves, closes and reopens the sandbox project, which releases AEDT's handles on the result
    files of deleted designs (a design re-created under a deleted design's name otherwise meets
    a modal 'file ... is in use' dialog on save); returns a fresh app for `design`.
    """
    desktop.odesktop.SetActiveProject(PROJECT_NAME).Save()
    desktop.odesktop.CloseProject(PROJECT_NAME)
    desktop.odesktop.OpenProject(str(PROJECT_PATH))
    app = Hfss(project=PROJECT_NAME, design=design, new_desktop=False, close_on_exit=False)
    check_sandbox(app)
    return app


def release(desktop):
    """Releases the session connection; AEDT and its projects stay open."""
    desktop.release_desktop(close_projects=False, close_on_exit=False)
