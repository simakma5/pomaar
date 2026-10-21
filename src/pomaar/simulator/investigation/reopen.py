#!/usr/bin/env python3
"""Saves, closes and reopens the investigation project in the GUI session (releases stale files)."""

from pomaar.simulator.investigation.project import open_project, release, reopen_project

if __name__ == "__main__":
    desktop, app = open_project("BistaticityAnalysis")
    try:
        app = reopen_project(desktop, app, "BistaticityAnalysis")
        print("Reopened; designs:", list(app.oproject.GetTopDesignList()))
    finally:
        release(desktop)
