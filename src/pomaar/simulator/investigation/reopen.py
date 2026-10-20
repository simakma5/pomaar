#!/usr/bin/env python3
"""Saves, closes and reopens the sandbox project in the GUI session (releases stale files)."""

from pomaar.simulator.investigation.sandbox import open_sandbox, release, reopen_sandbox

if __name__ == "__main__":
    desktop, _ = open_sandbox("BistaticityAnalysis")
    try:
        app = reopen_sandbox(desktop, "MultiRxTrihedral")
        print("Reopened; designs:", list(app.oproject.GetTopDesignList()))
    finally:
        release(desktop)
