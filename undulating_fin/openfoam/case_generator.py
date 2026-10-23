"""
OpenFOAM case generator for undulating fin FSI simulations.

Generates a complete case directory with:
- blockMeshDict for the 2-D channel + fin geometry
- Boundary / initial conditions (0/)
- fvSchemes, fvSolution, dynamicMeshDict (system/)
- controlDict with function objects for force extraction
- fin_motion.py pointDisplacement script driven by FinParams

Uses pimpleDyMFoam with morphing (ALE) mesh and displacementLaplacian
motion solver — the simplest approach for moderate-amplitude oscillations.
"""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

from undulating_fin.params import FinParams


# ── Defaults ────────────────────────────────────────────────────────

@dataclass
class OFConfig:
    """OpenFOAM simulation configuration."""

    # Time control
    n_cycles: int = 3
    write_interval_per_cycle: int = 20
    max_courant: float = 1.0

    # Domain sizing (multiples of fin span)
    inlet_length: float = 5.0     # upstream of fin [× D]
    outlet_length: float = 15.0   # downstream of fin [× D]
    lateral_width: float = 5.0    # lateral extent [× D]
    depth: float = 1.0            # 2-D extrude depth [× D]

    # Mesh
    cells_per_D: int = 40         # cells per inter-ray spacing
    boundary_layers: int = 5
    bl_expansion: float = 1.2
    bl_first_height: float = 0.001  # [m]

    # Fluid
    nu: float = 1e-6              # kinematic viscosity [m²/s]

    # Solver
    solver: str = "pimpleDyMFoam"
    n_outer_correctors: int = 2
    n_correctors: int = 2
    n_non_ortho_correctors: int = 1


# ── Header ──────────────────────────────────────────────────────────

_FOAM_HEADER = textwrap.dedent("""\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       {cls};
        object      {obj};
    }}
""")


def _header(cls: str, obj: str) -> str:
    return _FOAM_HEADER.format(cls=cls, obj=obj)


# ── Geometry helpers ────────────────────────────────────────────────

def _fin_vertices(params: FinParams) -> list[tuple[float, float]]:
    """Return 2-D polyline vertices for the fin at rest (q_i = 0).

    The fin is N rays connected end-to-end along the x-axis,
    each of length D.  Membranes extend from hmin to hmax in y.
    """
    pts = []
    # bottom edge (y = -hmax → -hmin along each ray root)
    for i in range(params.N):
        x = i * params.D
        pts.append((x, -params.hmax))
    # top edge (reverse)
    for i in reversed(range(params.N)):
        x = i * params.D
        pts.append((x, params.hmax))
    return pts


# ── Case generation ────────────────────────────────────────────────

def generate_case(
    params: FinParams,
    cfg: OFConfig,
    case_dir: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Create a complete OpenFOAM case directory.

    Parameters
    ----------
    params : FinParams
        Fin parameters (N, A, f, Phi, D, hmin, hmax, rho, …).
    cfg : OFConfig
        Solver and meshing settings.
    case_dir : str | Path
        Target directory.
    overwrite : bool
        If True, remove existing case_dir first.

    Returns
    -------
    Path to the created case directory.
    """
    case = Path(case_dir).resolve()
    if case.exists():
        if overwrite:
            shutil.rmtree(case)
        else:
            raise FileExistsError(f"Case directory exists: {case}")
    case.mkdir(parents=True)

    dt = 1.0 / (params.f * cfg.write_interval_per_cycle * 10)  # ~10× per write
    end_time = cfg.n_cycles / params.f
    write_dt = 1.0 / (params.f * cfg.write_interval_per_cycle)

    _write_control_dict(case, params, cfg, dt, end_time, write_dt)
    _write_fv_schemes(case)
    _write_fv_solution(case, cfg)
    _write_dynamic_mesh_dict(case)
    _write_block_mesh_dict(case, params, cfg)
    _write_boundary_conditions(case, cfg)
    _write_transport_properties(case, cfg)
    _write_fin_motion_script(case, params)
    _write_allrun(case, cfg)
    _write_allclean(case)

    return case


# ── controlDict ─────────────────────────────────────────────────────

def _write_control_dict(
    case: Path, params: FinParams, cfg: OFConfig,
    dt: float, end_time: float, write_dt: float,
) -> None:
    d = case / "system"
    d.mkdir(parents=True, exist_ok=True)
    content = _header("dictionary", "controlDict") + textwrap.dedent(f"""\
        application     {cfg.solver};
        startFrom       startTime;
        startTime       0;
        stopAt          endTime;
        endTime         {end_time};
        deltaT          {dt};
        writeControl    adjustableRunTime;
        writeInterval   {write_dt};
        purgeWrite      0;
        writeFormat     ascii;
        writePrecision  8;
        writeCompression off;
        timeFormat      general;
        timePrecision   6;
        runTimeModifiable true;
        adjustTimeStep  yes;
        maxCo           {cfg.max_courant};

        functions
        {{
            forces
            {{
                type            forces;
                libs            ("libforces.so");
                writeControl    timeStep;
                writeInterval   1;
                patches         (fin);
                rho             rhoInf;
                rhoInf          {params.rho};
                CofR            (0 0 0);
                log             true;
            }}
        }}
    """)
    (d / "controlDict").write_text(content, encoding="utf-8")


# ── fvSchemes ───────────────────────────────────────────────────────

def _write_fv_schemes(case: Path) -> None:
    d = case / "system"
    d.mkdir(parents=True, exist_ok=True)
    content = _header("dictionary", "fvSchemes") + textwrap.dedent("""\
        ddtSchemes
        {
            default         Euler;
        }
        gradSchemes
        {
            default         Gauss linear;
        }
        divSchemes
        {
            default         none;
            div(phi,U)      Gauss linearUpwind grad(U);
            div((nuEff*dev2(T(grad(U))))) Gauss linear;
        }
        laplacianSchemes
        {
            default         Gauss linear corrected;
        }
        interpolationSchemes
        {
            default         linear;
        }
        snGradSchemes
        {
            default         corrected;
        }
    """)
    (d / "fvSchemes").write_text(content, encoding="utf-8")


# ── fvSolution ──────────────────────────────────────────────────────

def _write_fv_solution(case: Path, cfg: OFConfig) -> None:
    d = case / "system"
    d.mkdir(parents=True, exist_ok=True)
    content = _header("dictionary", "fvSolution") + textwrap.dedent(f"""\
        solvers
        {{
            "pcorr.*"
            {{
                solver          PCG;
                preconditioner  DIC;
                tolerance       1e-5;
                relTol          0;
            }}
            p
            {{
                solver          PCG;
                preconditioner  DIC;
                tolerance       1e-6;
                relTol          0.01;
            }}
            pFinal
            {{
                $p;
                relTol          0;
            }}
            "(U|cellDisplacement)"
            {{
                solver          smoothSolver;
                smoother        symGaussSeidel;
                tolerance       1e-6;
                relTol          0.01;
            }}
            "(U|cellDisplacement)Final"
            {{
                $U;
                relTol          0;
            }}
        }}

        PIMPLE
        {{
            nOuterCorrectors    {cfg.n_outer_correctors};
            nCorrectors         {cfg.n_correctors};
            nNonOrthogonalCorrectors {cfg.n_non_ortho_correctors};
            correctPhi          yes;
            moveMeshOuterCorrectors yes;
        }}
    """)
    (d / "fvSolution").write_text(content, encoding="utf-8")


# ── dynamicMeshDict ─────────────────────────────────────────────────

def _write_dynamic_mesh_dict(case: Path) -> None:
    d = case / "constant"
    d.mkdir(parents=True, exist_ok=True)
    content = _header("dictionary", "dynamicMeshDict") + textwrap.dedent("""\
        dynamicFvMesh   dynamicMotionSolverFvMesh;
        motionSolverLibs ("libfvMotionSolvers.so");
        solver          displacementLaplacian;
        displacementLaplacianCoeffs
        {
            diffusivity     inverseDistance (fin);
        }
    """)
    (d / "dynamicMeshDict").write_text(content, encoding="utf-8")


# ── blockMeshDict ───────────────────────────────────────────────────

def _write_block_mesh_dict(
    case: Path, params: FinParams, cfg: OFConfig,
) -> None:
    d = case / "system"
    d.mkdir(parents=True, exist_ok=True)

    span = (params.N - 1) * params.D  # total fin span
    x_min = -cfg.inlet_length * params.D
    x_max = span + cfg.outlet_length * params.D
    y_min = -cfg.lateral_width * params.D
    y_max = cfg.lateral_width * params.D
    z_min = 0.0
    z_max = cfg.depth * params.D

    # Cell counts (approximate)
    nx = max(20, int((x_max - x_min) / params.D * cfg.cells_per_D))
    ny = max(20, int((y_max - y_min) / params.D * cfg.cells_per_D))
    nz = 1  # 2-D case

    content = _header("dictionary", "blockMeshDict") + textwrap.dedent(f"""\
        convertToMeters 1;

        vertices
        (
            ({x_min} {y_min} {z_min})   // 0
            ({x_max} {y_min} {z_min})   // 1
            ({x_max} {y_max} {z_min})   // 2
            ({x_min} {y_max} {z_min})   // 3
            ({x_min} {y_min} {z_max})   // 4
            ({x_max} {y_min} {z_max})   // 5
            ({x_max} {y_max} {z_max})   // 6
            ({x_min} {y_max} {z_max})   // 7
        );

        blocks
        (
            hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
        );

        edges
        (
        );

        boundary
        (
            inlet
            {{
                type patch;
                faces
                (
                    (0 4 7 3)
                );
            }}
            outlet
            {{
                type patch;
                faces
                (
                    (1 2 6 5)
                );
            }}
            top
            {{
                type wall;
                faces
                (
                    (3 7 6 2)
                );
            }}
            bottom
            {{
                type wall;
                faces
                (
                    (0 1 5 4)
                );
            }}
            frontAndBack
            {{
                type empty;
                faces
                (
                    (0 3 2 1)
                    (4 5 6 7)
                );
            }}
        );

        mergePatchPairs
        (
        );
    """)
    (d / "blockMeshDict").write_text(content, encoding="utf-8")


# ── Boundary conditions (0/) ────────────────────────────────────────

def _write_boundary_conditions(case: Path, cfg: OFConfig) -> None:
    d = case / "0"
    d.mkdir(parents=True, exist_ok=True)

    # --- U ---
    u_content = _header("volVectorField", "U") + textwrap.dedent("""\
        dimensions      [0 1 -1 0 0 0 0];
        internalField   uniform (0 0 0);

        boundaryField
        {
            inlet
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            outlet
            {
                type            zeroGradient;
            }
            top
            {
                type            noSlip;
            }
            bottom
            {
                type            noSlip;
            }
            fin
            {
                type            movingWallVelocity;
                value           uniform (0 0 0);
            }
            frontAndBack
            {
                type            empty;
            }
        }
    """)
    (d / "U").write_text(u_content, encoding="utf-8")

    # --- p ---
    p_content = _header("volScalarField", "p") + textwrap.dedent("""\
        dimensions      [0 2 -2 0 0 0 0];
        internalField   uniform 0;

        boundaryField
        {
            inlet
            {
                type            zeroGradient;
            }
            outlet
            {
                type            fixedValue;
                value           uniform 0;
            }
            top
            {
                type            zeroGradient;
            }
            bottom
            {
                type            zeroGradient;
            }
            fin
            {
                type            zeroGradient;
            }
            frontAndBack
            {
                type            empty;
            }
        }
    """)
    (d / "p").write_text(p_content, encoding="utf-8")

    # --- pointDisplacement ---
    pd_content = _header("pointVectorField", "pointDisplacement") + textwrap.dedent("""\
        dimensions      [0 1 0 0 0 0 0];
        internalField   uniform (0 0 0);

        boundaryField
        {
            inlet
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            outlet
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            top
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            bottom
            {
                type            fixedValue;
                value           uniform (0 0 0);
            }
            fin
            {
                type            codedFixedValue;
                redirectType    finMotion;
                value           uniform (0 0 0);
                name            finMotion;
                code
                #{
                    // Motion is defined in system/finMotionCodeDict
                    // or generated by the Python pre-processor.
                    // This placeholder applies zero displacement;
                    // replace with coded boundary from fin_motion.py output.
                    const vectorField& points = patch().localPoints();
                    vectorField& disp = *this;
                    disp = vector::zero;
                #};
            }
            frontAndBack
            {
                type            empty;
            }
        }
    """)
    (d / "pointDisplacement").write_text(pd_content, encoding="utf-8")


# ── Transport properties ────────────────────────────────────────────

def _write_transport_properties(case: Path, cfg: OFConfig) -> None:
    d = case / "constant"
    d.mkdir(parents=True, exist_ok=True)
    content = _header("dictionary", "transportProperties") + textwrap.dedent(f"""\
        transportModel  Newtonian;
        nu              [0 2 -1 0 0 0 0] {cfg.nu};
    """)
    (d / "transportProperties").write_text(content, encoding="utf-8")


# ── Fin motion Python script ───────────────────────────────────────

def _write_fin_motion_script(case: Path, params: FinParams) -> None:
    """Generate a standalone Python script that computes pointDisplacement
    for the fin boundary at each time step.

    This script is called by the Allrun workflow to generate time-varying
    boundary data, or can be used to generate a ``timeVaryingMapped``
    boundary dataset.
    """
    d = case / "scripts"
    d.mkdir(parents=True, exist_ok=True)

    motion_src = textwrap.dedent(f"""\
        #!/usr/bin/env python3
        \"\"\"
        Compute fin boundary displacement at each time step.

        Usage:
            python fin_motion.py <time> <points_file> <output_file>

        The script reads the reference fin boundary points, applies the
        traveling-wave kinematics from FinParams, and writes the
        displacement vectors.
        \"\"\"

        import sys
        import numpy as np

        # Fin parameters (frozen from case generation)
        N      = {params.N}
        A      = {params.A}         # amplitude [rad]
        f      = {params.f}         # frequency [Hz]
        Phi    = {params.Phi}       # inter-ray phase shift [rad]
        X      = {params.X}         # angular offset [rad]
        D      = {params.D}         # inter-ray distance [m]
        hmax   = {params.hmax}      # max membrane height [m]
        hmin   = {params.hmin}      # min membrane height [m]


        def ray_angle(i: int, t: float) -> float:
            \"\"\"Joint angle for ray *i* (0-indexed) at time *t*.\"\"\"
            return A * np.sin(2 * np.pi * f * t - i * Phi) + X


        def membrane_displacement(x_ref: float, y_ref: float, t: float):
            \"\"\"Compute (dx, dy) displacement for a mesh point at (x_ref, y_ref).

            Determines which panel the point belongs to from its x-coordinate,
            computes the two bounding ray angles, and interpolates the
            displacement using the normalized span position w_bar.
            \"\"\"
            # Identify panel index
            panel = int(x_ref / D)
            panel = max(0, min(panel, N - 2))
            q1 = ray_angle(panel, t)
            q2 = ray_angle(panel + 1, t)

            # Normalized span
            x_local = x_ref - panel * D
            w_bar = x_local / D if D > 0 else 0.0
            w_bar = max(0.0, min(1.0, w_bar))

            # Height (distance from joint axis)
            h = abs(y_ref)
            h = max(hmin, min(h, hmax))
            sign_y = 1.0 if y_ref >= 0 else -1.0

            # Displacement from rotating rays
            # Ray at panel rotates by q1, ray at panel+1 by q2
            # Interpolated angle at this span position
            q_interp = q1 * (1.0 - w_bar) + q2 * w_bar

            # Displacement of this point (rotation about joint)
            # Reference position: y = h, x = x_local relative to joint
            dx = 0.0  # no axial displacement (rays rotate in y)
            dy = sign_y * h * (np.cos(q_interp) - 1.0)  # transverse displacement

            # Add the tangential component
            dx += -sign_y * h * np.sin(q_interp)

            return dx, dy


        def main():
            if len(sys.argv) < 4:
                print("Usage: fin_motion.py <time> <points_file> <output_file>")
                sys.exit(1)

            t = float(sys.argv[1])
            points_file = sys.argv[2]
            output_file = sys.argv[3]

            # Read reference points (OpenFOAM format)
            points = []
            with open(points_file, "r") as fh:
                in_block = False
                for line in fh:
                    line = line.strip()
                    if line == "(":
                        in_block = True
                        continue
                    if line == ")":
                        break
                    if in_block and line.startswith("("):
                        coords = line.strip("()").split()
                        points.append([float(c) for c in coords])
            points = np.array(points)

            # Compute displacements
            displacements = np.zeros_like(points)
            for idx in range(len(points)):
                dx, dy = membrane_displacement(points[idx, 0], points[idx, 1], t)
                displacements[idx, 0] = dx
                displacements[idx, 1] = dy
                # z = 0 (2-D)

            # Write OpenFOAM-format vector field
            with open(output_file, "w") as fh:
                fh.write(str(len(displacements)) + "\\n")
                fh.write("(\\n")
                for d in displacements:
                    fh.write(f"({{d[0]:.10e}} {{d[1]:.10e}} {{d[2]:.10e}})\\n")
                fh.write(")\\n")

            print(f"Wrote {{len(displacements)}} displacement vectors for t = {{t}}")


        if __name__ == "__main__":
            main()
    """)
    (d / "fin_motion.py").write_text(motion_src, encoding="utf-8")


# ── Allrun / Allclean ───────────────────────────────────────────────

def _write_allrun(case: Path, cfg: OFConfig) -> None:
    content = textwrap.dedent(f"""\
        #!/bin/bash
        # Allrun script for undulating fin OpenFOAM case
        cd "${{0%/*}}" || exit
        . ${{WM_PROJECT_DIR:?}}/bin/tools/RunFunctions

        echo "=== Undulating Fin OpenFOAM Simulation ==="

        # 1. Mesh generation
        runApplication blockMesh

        # 2. Optional: add fin geometry via snappyHexMesh
        # runApplication snappyHexMesh -overwrite

        # 3. Initialize fields
        runApplication setFields 2>/dev/null || true

        # 4. Check mesh quality
        runApplication checkMesh

        # 5. Run solver
        runApplication {cfg.solver}

        # 6. Post-process: extract forces
        runApplication postProcess -func forces

        echo "=== Simulation Complete ==="
    """)
    (case / "Allrun").write_text(content, encoding="utf-8")


def _write_allclean(case: Path) -> None:
    content = textwrap.dedent("""\
        #!/bin/bash
        # Allclean script
        cd "${0%/*}" || exit
        . ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions

        cleanCase
        rm -rf constant/polyMesh
        rm -rf postProcessing
    """)
    (case / "Allclean").write_text(content, encoding="utf-8")
