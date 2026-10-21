#!/usr/bin/env python3
"""
Lists the projects open in the live GUI session, with their folders and designs (read-only).

    ansys-vnc exec python -m pomaar.agent.list_projects
"""

from pomaar.simulator.hfss_array_builder import connect_live_session

if __name__ == "__main__":
    desktop = connect_live_session()
    try:
        for name in desktop.project_list:
            project = desktop.odesktop.SetActiveProject(name)
            print(name, project.GetPath(), list(project.GetTopDesignList()))
    finally:
        desktop.release_desktop(close_projects=False, close_on_exit=False)
