#!/usr/bin/env python3
"""
Sets the initial (SBR+ surface) mesh settings of a sandbox design: the surface-approximation
slider (1 coarse ... 10 fine) and curvilinear elements, which decide how finely curved targets
such as the sphere are faceted:

    ansys-vnc exec python src/pomaar/simulator/investigation/mesh.py MultiRxSphere \
        --slider 10 [--curvilinear]
"""

import argparse

from pomaar.simulator.investigation.sandbox import open_sandbox, release


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("design")
    parser.add_argument("--slider", type=int, default=5, help="1 (coarse) ... 10 (fine)")
    parser.add_argument("--curvilinear", action="store_true")
    args = parser.parse_args(argv)

    desktop, app = open_sandbox(args.design)
    try:
        app.odesign.GetModule("MeshSetup").InitialMeshSettings(
            [
                "NAME:MeshSettings",
                [
                    "NAME:GlobalSurfApproximation",
                    "CurvedSurfaceApproxChoice:=", "UseSlider",
                    "SliderMeshSettings:=", args.slider,
                ],
                ["NAME:GlobalCurvilinear", "Apply:=", args.curvilinear],
                ["NAME:GlobalModelRes", "UseAutoLength:=", True],
                "MeshMethod:=", "Auto",
                "UseTolerantSBRMesh:=", True,
                "UseDynamicMeshing:=", False,
                "DynamicSurfaceResolution:=", False,
                "UseFlexMeshingForTAUvolumeMesh:=", False,
                "UseAlternativeMeshMethodsAsFallBack:=", True,
                "AllowPhiForLayeredGeometry:=", True,
            ]
        )  # fmt: skip
        app.save_project()
        print(f"{args.design}: slider {args.slider}, curvilinear {args.curvilinear}")
    finally:
        release(desktop)


if __name__ == "__main__":
    main()
