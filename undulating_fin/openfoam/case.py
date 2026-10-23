r"""
OpenFOAM case generation for the undulating fin.

Builds a complete, runnable 3-D case: background block mesh,
``snappyHexMesh`` around the fin surface, prescribed membrane motion on
a morphing mesh, and force extraction on the fin patch.

Every dictionary here follows a pattern that was executed successfully
under OpenFOAM v2412 before being written — see ``reference/README.md``
for the transcripts.  The predecessor of this module emitted cases that
had never been run, and all of them were broken: the mesh contained no
fin at all, the motion boundary condition was a stub that applied zero
displacement, and the solver it named had been removed from OpenFOAM
years earlier.

The case is 3-D because the reference COMSOL simulation is: the fin's
rays stack along *z* and rotate in the *x*–*y* plane, so a 2-D domain
cannot produce three independent force components.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from undulating_fin.params import FinParams

from .geometry import SPAN_SIGN
from .stl import build_fin_solid, write_stl_ascii

__all__ = ["CaseConfig", "generate_case"]


# ── Configuration ───────────────────────────────────────────────────

@dataclass
class CaseConfig:
    """Numerical setup, independent of the fin's physical parameters."""

    # ── Time ──
    n_cycles: float = 3.0
    write_per_cycle: int = 40
    dt: float = 5.0e-4
    adjust_time_step: bool = True
    max_courant: float = 5.0

    # ── Domain, as multiples of the fin's radial reach (hmax) ──
    radial_extent: float = 4.0
    span_upstream: float = 3.0     # ×D beyond the first ray
    span_downstream: float = 5.0   # ×D beyond the last ray

    # ── Mesh ──
    background_cell: float = 0.04  # [m]
    surface_refinement: int = 4    # snappy levels on the fin
    volume_refinement: int = 2     # levels in the box around the fin
    n_layers: int = 0              # prism layers (0 = none; see note)
    membrane_thickness: float = 6.0e-3

    # ── Fluid ──
    rho: float = 1000.0            # [kg/m³] — matches the COMSOL run
    nu: float = 1.0e-6             # [m²/s]

    # ── Solver ──
    n_outer_correctors: int = 2
    n_correctors: int = 2
    n_non_orthogonal_correctors: int = 1
    # 4 ranks measured fastest on this host; 6 is slower than 4 and
    # `mpirun -np 8` fails, as OpenMPI counts physical cores (7), not
    # threads (14).
    n_procs: int = 4

    # ── Surface tessellation ──
    stl_n_h: int = 60
    stl_n_w: int = 24

    def __post_init__(self) -> None:
        if self.membrane_thickness <= 0:
            raise ValueError("membrane_thickness must be positive")
        if self.n_procs < 1:
            raise ValueError("n_procs must be at least 1")
        finest = self.background_cell / 2 ** self.surface_refinement
        if finest > 0.5 * self.membrane_thickness:
            raise ValueError(
                f"finest cell {finest*1e3:.2f} mm cannot resolve a "
                f"{self.membrane_thickness*1e3:.2f} mm membrane; raise "
                "surface_refinement, shrink background_cell, or thicken "
                "the membrane"
            )


_HEADER = textwrap.dedent(
    """\
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       {cls};
        object      {obj};
    }}

    """
)


def _write(path: Path, cls: str, obj: str, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = _HEADER.format(cls=cls, obj=obj) + textwrap.dedent(body).lstrip("\n")
    path.write_text(text, encoding="utf-8", newline="\n")


# ── Domain ──────────────────────────────────────────────────────────

def _domain_bounds(p: FinParams, cfg: CaseConfig) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned domain box.

    The fin sweeps a fan of radius ``hmax`` in the *x*–*y* plane and
    occupies ``z ∈ [-(N-1)D, 0]`` because the rays stack along −z.
    """
    r = cfg.radial_extent * p.hmax
    span = (p.N - 1) * p.D

    lo = np.array([-r, -r, -span - cfg.span_downstream * p.D])
    hi = np.array([+r, +r, cfg.span_upstream * p.D])
    return lo, hi


def _inside_point(p: FinParams, cfg: CaseConfig) -> np.ndarray:
    """A point in the fluid, well clear of the fin.

    snappyHexMesh keeps the region connected to this point.  Placing it
    near a domain corner keeps it outside the swept fan for any phase.
    """
    lo, hi = _domain_bounds(p, cfg)
    return lo + 0.05 * (hi - lo)


# ── Dictionaries ────────────────────────────────────────────────────

def _block_mesh(case: Path, p: FinParams, cfg: CaseConfig) -> tuple[int, int, int]:
    lo, hi = _domain_bounds(p, cfg)
    n = np.maximum(np.round((hi - lo) / cfg.background_cell).astype(int), 1)
    (x0, y0, z0), (x1, y1, z1) = lo, hi

    _write(
        case / "system/blockMeshDict", "dictionary", "blockMeshDict",
        f"""
        scale   1;

        vertices
        (
            ({x0:.6f} {y0:.6f} {z0:.6f})
            ({x1:.6f} {y0:.6f} {z0:.6f})
            ({x1:.6f} {y1:.6f} {z0:.6f})
            ({x0:.6f} {y1:.6f} {z0:.6f})
            ({x0:.6f} {y0:.6f} {z1:.6f})
            ({x1:.6f} {y0:.6f} {z1:.6f})
            ({x1:.6f} {y1:.6f} {z1:.6f})
            ({x0:.6f} {y1:.6f} {z1:.6f})
        );

        blocks
        (
            hex (0 1 2 3 4 5 6 7) ({n[0]} {n[1]} {n[2]}) simpleGrading (1 1 1)
        );

        edges ();

        boundary
        (
            farfield
            {{
                type patch;
                faces
                (
                    (0 4 7 3) (1 2 6 5) (0 1 5 4)
                    (3 7 6 2) (0 3 2 1) (4 5 6 7)
                );
            }}
        );

        mergePatchPairs ();
        """,
    )
    return tuple(int(v) for v in n)


def _snappy(case: Path, p: FinParams, cfg: CaseConfig) -> None:
    lo, hi = _domain_bounds(p, cfg)
    ip = _inside_point(p, cfg)

    # Refinement box: the fan the fin sweeps, with a margin.
    margin = 0.25 * p.hmax
    span = (p.N - 1) * p.D
    bl = np.array([-p.hmax - margin, -p.hmax - margin, -span - margin])
    bh = np.array([+p.hmax + margin, +p.hmax + margin, +margin])

    layers = (
        ""
        if cfg.n_layers <= 0
        else f"""
            fin
            {{
                nSurfaceLayers {cfg.n_layers};
            }}
        """
    )

    _write(
        case / "system/snappyHexMeshDict", "dictionary", "snappyHexMeshDict",
        f"""
        castellatedMesh true;
        snap            true;
        addLayers       {"true" if cfg.n_layers > 0 else "false"};

        geometry
        {{
            fin.stl
            {{
                type triSurfaceMesh;
                name fin;
            }}
            refinementBox
            {{
                type searchableBox;
                min ({bl[0]:.6f} {bl[1]:.6f} {bl[2]:.6f});
                max ({bh[0]:.6f} {bh[1]:.6f} {bh[2]:.6f});
            }}
        }}

        castellatedMeshControls
        {{
            maxLocalCells       2000000;
            maxGlobalCells      8000000;
            minRefinementCells  10;
            maxLoadUnbalance    0.10;
            nCellsBetweenLevels 3;

            features
            (
                {{ file "fin.eMesh"; level {cfg.surface_refinement}; }}
            );

            refinementSurfaces
            {{
                fin
                {{
                    level ({cfg.surface_refinement} {cfg.surface_refinement});
                    patchInfo {{ type wall; }}
                }}
            }}

            resolveFeatureAngle 30;

            refinementRegions
            {{
                refinementBox
                {{
                    mode inside;
                    levels ((1e15 {cfg.volume_refinement}));
                }}
            }}

            locationInMesh ({ip[0]:.6f} {ip[1]:.6f} {ip[2]:.6f});
            allowFreeStandingZoneFaces false;
        }}

        snapControls
        {{
            nSmoothPatch    3;
            tolerance       2.0;
            nSolveIter      50;
            nRelaxIter      5;
            nFeatureSnapIter 10;
            implicitFeatureSnap    false;
            explicitFeatureSnap    true;
            multiRegionFeatureSnap false;
        }}

        addLayersControls
        {{
            relativeSizes    true;
            expansionRatio   1.2;
            finalLayerThickness 0.4;
            minThickness     0.05;
            nGrow            0;
            featureAngle     120;
            nRelaxIter       5;
            nSmoothSurfaceNormals 1;
            nSmoothNormals   3;
            nSmoothThickness 10;
            maxFaceThicknessRatio 0.5;
            maxThicknessToMedialRatio 0.3;
            minMedialAxisAngle 90;
            nBufferCellsNoExtrude 0;
            nLayerIter       50;
            layers
            {{{layers}}}
        }}

        meshQualityControls
        {{
            maxNonOrtho         65;
            maxBoundarySkewness 20;
            maxInternalSkewness 4;
            maxConcave          80;
            minVol              1e-13;
            minTetQuality       1e-15;
            minArea             -1;
            minTwist            0.02;
            minDeterminant      0.001;
            minFaceWeight       0.05;
            minVolRatio         0.01;
            minTriangleTwist    -1;
            nSmoothScale        4;
            errorReduction      0.75;
        }}

        mergeTolerance 1e-6;
        """,
    )

    _write(
        case / "system/surfaceFeatureExtractDict", "dictionary",
        "surfaceFeatureExtractDict",
        """
        fin.stl
        {
            extractionMethod    extractFromSurface;
            includedAngle       150;
            subsetFeatures      { nonManifoldEdges no; openEdges yes; }
            writeObj            no;
        }
        """,
    )


def _control(case: Path, p: FinParams, cfg: CaseConfig) -> None:
    end_time = cfg.n_cycles / p.f
    write_interval = 1.0 / (p.f * cfg.write_per_cycle)

    _write(
        case / "system/controlDict", "dictionary", "controlDict",
        f"""
        application     pimpleFoam;

        startFrom       startTime;
        startTime       0;
        stopAt          endTime;
        endTime         {end_time:.6f};

        deltaT          {cfg.dt:.8g};
        adjustTimeStep  {"yes" if cfg.adjust_time_step else "no"};
        maxCo           {cfg.max_courant};
        maxDeltaT       {cfg.dt:.8g};

        writeControl    adjustableRunTime;
        writeInterval   {write_interval:.6f};
        purgeWrite      0;
        writeFormat     ascii;
        writePrecision  8;
        writeCompression on;
        timeFormat      general;
        timePrecision   8;
        runTimeModifiable false;

        functions
        {{
            forces
            {{
                type            forces;
                libs            (forces);
                writeControl    timeStep;
                writeInterval   1;
                patches         (fin);
                rho             rhoInf;
                rhoInf          {cfg.rho};
                CofR            (0 0 0);
                writeFields     no;
                log             no;
            }}

            meshQuality
            {{
                // Morphing meshes degrade as the fin swings; this makes
                // the degradation visible instead of letting it end the
                // run with an opaque solver failure.
                type            cellQuality;
                libs            (fieldFunctionObjects);
                writeControl    writeTime;
            }}

            minMaxU
            {{
                type            fieldMinMax;
                libs            (fieldFunctionObjects);
                writeControl    timeStep;
                writeInterval   20;
                fields          (U p);
                log             yes;
            }}
        }}
        """,
    )


def _fv_schemes(case: Path) -> None:
    _write(
        case / "system/fvSchemes", "dictionary", "fvSchemes",
        """
        ddtSchemes
        {
            // Second order, bounded — the flow is oscillatory and a
            // first-order scheme would damp the very forces we measure.
            default         backward;
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
            default         Gauss linear limited corrected 0.33;
        }

        interpolationSchemes
        {
            default         linear;
        }

        snGradSchemes
        {
            default         limited corrected 0.33;
        }

        wallDist
        {
            method          meshWave;
        }
        """,
    )


def _fv_solution(case: Path, cfg: CaseConfig) -> None:
    # cellDisplacement is a symmetric (Laplacian) matrix and must not be
    # folded into the U regex — DILU is rejected for symmetric matrices.
    _write(
        case / "system/fvSolution", "dictionary", "fvSolution",
        f"""
        solvers
        {{
            "(p|pcorr).*"
            {{
                solver          GAMG;
                tolerance       1e-7;
                relTol          0.01;
                smoother        GaussSeidel;
            }}

            pFinal
            {{
                $p;
                relTol          0;
            }}

            "(U|UFinal)"
            {{
                solver          PBiCGStab;
                preconditioner  DILU;
                tolerance       1e-8;
                relTol          0;
            }}

            "cellDisplacement.*"
            {{
                solver          PCG;
                preconditioner  DIC;
                tolerance       1e-8;
                relTol          0;
            }}
        }}

        PIMPLE
        {{
            correctPhi                  yes;
            nOuterCorrectors            {cfg.n_outer_correctors};
            nCorrectors                 {cfg.n_correctors};
            nNonOrthogonalCorrectors    {cfg.n_non_orthogonal_correctors};
        }}

        relaxationFactors
        {{
            equations
            {{
                ".*"    1;
            }}
        }}
        """,
    )


def _constant(case: Path, cfg: CaseConfig) -> None:
    _write(
        case / "constant/transportProperties", "dictionary",
        "transportProperties",
        f"""
        transportModel  Newtonian;
        nu              {cfg.nu:.6g};
        """,
    )
    # ESI OpenFOAM expects turbulenceProperties; momentumTransport is
    # the OpenFOAM Foundation spelling and is not read here.
    _write(
        case / "constant/turbulenceProperties", "dictionary",
        "turbulenceProperties",
        """
        simulationType  laminar;
        """,
    )
    # motionSolverLibs is required: without it displacementLaplacian is
    # reported as an unknown solver type.
    _write(
        case / "constant/dynamicMeshDict", "dictionary", "dynamicMeshDict",
        """
        dynamicFvMesh       dynamicMotionSolverFvMesh;
        motionSolverLibs    (fvMotionSolvers);
        motionSolver        displacementLaplacian;

        displacementLaplacianCoeffs
        {
            // Diffusivity falls away from the fin, so cells near the
            // membrane move almost rigidly with it and the distortion
            // is pushed out into the coarse far field.
            diffusivity     inverseDistance (fin);
        }
        """,
    )


def _motion_kernel(p: FinParams) -> str:
    """Self-contained C++ for the membrane displacement.

    Emitted into ``codeInclude`` (file scope) rather than inline in the
    boundary condition body, so the arithmetic can be extracted between
    the markers, compiled on its own and checked against the Python
    reference — see ``test_cpp_transcription_matches_python``. A
    transcription that has silently drifted from its reference is
    exactly how the previous generator ended up with dead code.

    A direct transcription of :func:`geometry.membrane_displacement`.
    """
    return f"""// --- fin motion kernel begin ---
static inline void finMembraneDisp
(
    const double t,
    const double Px, const double Py, const double Pz,
    double& dx, double& dy
)
{{
    const int    N        = {p.N};
    const double A        = {p.A:.17g};
    const double freq     = {p.f:.17g};
    const double Phi      = {p.Phi:.17g};
    const double X        = {p.X:.17g};
    const double D        = {p.D:.17g};
    const double tRef     = 0.0;
    const double spanSign = {SPAN_SIGN:.1f};
    const double twoPi    = 6.283185307179586476925286766559;

    // Spanwise station: z = spanSign*(panel + w)*D, so the panel and
    // the normalised span coordinate both follow from z alone.
    double s = spanSign * Pz / D;
    int panel = (int)std::floor(s);
    if (panel < 0)     panel = 0;
    if (panel > N - 2) panel = N - 2;
    double w = s - (double)panel;
    if (w < 0.0) w = 0.0;
    if (w > 1.0) w = 1.0;

    const double q1r = A*std::sin(twoPi*freq*tRef - panel*Phi) + X;
    const double q2r = A*std::sin(twoPi*freq*tRef - (panel+1)*Phi) + X;
    const double q1  = A*std::sin(twoPi*freq*t    - panel*Phi) + X;
    const double q2  = A*std::sin(twoPi*freq*t    - (panel+1)*Phi) + X;

    // Blended radial direction; not a unit vector once the bounding
    // rays differ.
    const double uxr = (1.0-w)*std::cos(q1r) + w*std::cos(q2r);
    const double uyr = (1.0-w)*std::sin(q1r) + w*std::sin(q2r);
    const double uxt = (1.0-w)*std::cos(q1)  + w*std::cos(q2);
    const double uyt = (1.0-w)*std::sin(q1)  + w*std::sin(q2);

    double gRef = std::sqrt(uxr*uxr + uyr*uyr);
    if (gRef < 1e-300) gRef = 1e-300;
    const double gNow = std::sqrt(uxt*uxt + uyt*uyt);

    // Rotation by the swing of that direction plus a stretch by the
    // change in its length.  Both depend only on the station, so this
    // carries the inflated solid and not just the mid-surface.
    const double dphi  = std::atan2(uyt, uxt) - std::atan2(uyr, uxr);
    const double scale = gNow/gRef;
    const double cd = std::cos(dphi);
    const double sd = std::sin(dphi);

    const double x1 = scale*(Px*cd - Py*sd);
    const double y1 = scale*(Px*sd + Py*cd);

    dx = x1 - Px;
    dy = y1 - Py;
}}
// --- fin motion kernel end ---"""


def _fields(case: Path, p: FinParams, cfg: CaseConfig) -> None:
    _write(
        case / "0.orig/U", "volVectorField", "U",
        """
        dimensions      [0 1 -1 0 0 0 0];
        internalField   uniform (0 0 0);

        boundaryField
        {
            farfield
            {
                // Quiescent fluid: let the fin's own wake leave without
                // pinning the far field to zero velocity.
                type            pressureInletOutletVelocity;
                value           uniform (0 0 0);
            }

            fin
            {
                // Accounts for the flux through the moving face, which
                // a plain noSlip on a deforming mesh would not.
                type            movingWallVelocity;
                value           uniform (0 0 0);
            }
        }
        """,
    )

    _write(
        case / "0.orig/p", "volScalarField", "p",
        """
        dimensions      [0 2 -2 0 0 0 0];
        internalField   uniform 0;

        boundaryField
        {
            farfield
            {
                type            totalPressure;
                p0              uniform 0;
                value           uniform 0;
            }

            fin
            {
                type            zeroGradient;
            }
        }
        """,
    )

    kernel = _motion_kernel(p)
    _write(
        case / "0.orig/pointDisplacement", "pointVectorField", "pointDisplacement",
        f"""
        dimensions      [0 1 0 0 0 0 0];
        internalField   uniform (0 0 0);

        boundaryField
        {{
            farfield
            {{
                type            fixedValue;
                value           uniform (0 0 0);
            }}

            fin
            {{
                type            codedFixedValue;
                value           uniform (0 0 0);
                name            finMembraneMotion;

                codeInclude
                #{{
                    #include "points0MotionSolver.H"
                    #include <cmath>
{textwrap.indent(kernel, " " * 20)}
                #}};

                codeOptions
                #{{
                    -I$(LIB_SRC)/dynamicMesh/lnInclude
                #}};

                codeLibs
                #{{
                    -ldynamicMesh
                #}};

                code
                #{{
                    const pointPatch& pp = this->patch();
                    const polyMesh& m = pp.boundaryMesh().mesh()();
                    const points0MotionSolver& ms =
                        m.lookupObject<points0MotionSolver>("dynamicMeshDict");
                    const pointField& p0 = ms.points0();
                    const labelList& mp = pp.meshPoints();
                    const scalar t = this->db().time().value();

                    vectorField d(mp.size(), Zero);
                    forAll(mp, i)
                    {{
                        const point& P = p0[mp[i]];
                        double dx = 0.0, dy = 0.0;
                        finMembraneDisp(t, P.x(), P.y(), P.z(), dx, dy);
                        d[i] = vector(dx, dy, 0.0);
                    }}
                    operator==(d);
                #}};
            }}
        }}
        """,
    )


def _decompose(case: Path, cfg: CaseConfig) -> None:
    _write(
        case / "system/decomposeParDict", "dictionary", "decomposeParDict",
        f"""
        numberOfSubdomains  {cfg.n_procs};
        method              scotch;
        """,
    )


def _scripts(case: Path, cfg: CaseConfig) -> None:
    par = cfg.n_procs > 1
    run = f"""\
#!/bin/bash
cd "${{0%/*}}" || exit 1
. "${{WM_PROJECT_DIR:?set up OpenFOAM first}}/bin/tools/RunFunctions"

# Mesh
runApplication surfaceFeatureExtract
runApplication blockMesh
runApplication snappyHexMesh -overwrite
runApplication checkMesh -constant

# snappyHexMesh names the fin patch itself; the fields are staged only
# now so that they are mapped onto the final mesh.
restore0Dir
"""
    if par:
        run += f"""
runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructPar
"""
    else:
        run += """
runApplication $(getApplication)
"""
    (case / "Allrun").write_text(run, encoding="utf-8", newline="\n")
    (case / "Allrun").chmod(0o755)

    clean = """\
#!/bin/bash
cd "${0%/*}" || exit 1
. "${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions"

cleanCase0
rm -rf dynamicCode postProcessing constant/extendedFeatureEdgeMesh
rm -f constant/triSurface/fin.eMesh
"""
    (case / "Allclean").write_text(clean, encoding="utf-8", newline="\n")
    (case / "Allclean").chmod(0o755)


# ── Entry point ─────────────────────────────────────────────────────

def generate_case(
    params: FinParams,
    out_dir: str | Path,
    cfg: CaseConfig | None = None,
) -> Path:
    """Write a complete OpenFOAM case and return its directory.

    Parameters
    ----------
    params
        Fin geometry and kinematics.
    out_dir
        Directory to create.  Overwritten if it exists.
    cfg
        Numerical setup; defaults to :class:`CaseConfig`.
    """
    cfg = cfg or CaseConfig()
    case = Path(out_dir)
    if case.exists():
        import shutil

        shutil.rmtree(case)
    case.mkdir(parents=True)

    verts, faces = build_fin_solid(
        params, t=0.0, n_h=cfg.stl_n_h, n_w=cfg.stl_n_w,
        thickness=cfg.membrane_thickness,
    )
    write_stl_ascii(case / "constant/triSurface/fin.stl", verts, faces, name="fin")

    _block_mesh(case, params, cfg)
    _snappy(case, params, cfg)
    _control(case, params, cfg)
    _fv_schemes(case)
    _fv_solution(case, cfg)
    _constant(case, cfg)
    _fields(case, params, cfg)
    _decompose(case, cfg)
    _scripts(case, cfg)

    return case
