# Verified OpenFOAM patterns

Dictionaries in this directory were **run successfully** under OpenFOAM
v2412 before being saved. The case generator is built on them rather
than on patterns assumed to work — the previous generator produced
cases that had never been executed, and every one of them was broken.

Each `.verified` file is the exact text that ran, not a reconstruction.

## `pointDisplacement.verified` — prescribing membrane motion

The fin's motion is a function of each mesh point's **reference**
position, so the boundary condition needs `points0`, not the current
points (`patch().localPoints()` returns the latter and would feed back
on itself).

`points0` is registered with `IOobject::NO_REGISTER`, so it cannot be
looked up by name from the object registry. It has to be reached
through the motion solver, which *is* registered — as an `IOdictionary`
under the name `dynamicMeshDict`:

```cpp
const pointPatch& pp = this->patch();
const polyMesh& m = pp.boundaryMesh().mesh()();
const points0MotionSolver& ms =
    m.lookupObject<points0MotionSolver>("dynamicMeshDict");
const pointField& p0 = ms.points0();
const labelList& mp = pp.meshPoints();      // patch point -> mesh point
```

with

```
codeInclude  #{ #include "points0MotionSolver.H" #};
codeOptions  #{ -I$(LIB_SRC)/dynamicMesh/lnInclude #};
codeLibs     #{ -ldynamicMesh #};
```

**Verification.** The probe imposed `d = 0.05·y₀·sin(100t)` and reported
`p0.size() = 729` (9³ points of an 8³-cell block), `patchPts = 81` (9²),
and displacements matching the closed form at both time steps:

| t | expected | reported |
|---|---|---|
| 0.001 | 0.05·0.1·sin(0.1) = 4.99167e-4 | 4.99167e-4 |
| 0.002 | 0.05·0.1·sin(0.2) = 9.93347e-4 | 9.93347e-4 |

Confirming the values, not merely that it compiled, is what rules out
reading the wrong point set.

## `fvSolution.verified` — solver entries a moving mesh needs

Four entries that a static case does not require, each of which was hit
as a fatal error while getting the probe to run:

| Entry | Why |
|---|---|
| `cellDisplacement.*` | the mesh-motion Laplacian. It is **symmetric**: use `PCG`/`DIC`. Folding it into the `U` regex picks `PBiCGStab`/`DILU` and dies with *"Unknown symmetric matrix preconditioner type DILU"* |
| `pcorr.*` | PIMPLE's mesh-motion flux correction, solved only when the mesh moves |
| `*Final` variants | the regexes must cover `UFinal`, `pFinal`, `cellDisplacementFinal`, `pcorrFinal` |
| `pRefCell` / `pRefValue` | only for a fully closed domain; unnecessary once the case has an outlet |

## Other gotchas confirmed the same way

- `motionSolverLibs (fvMotionSolvers);` is required in
  `dynamicMeshDict`, otherwise `displacementLaplacian` is reported as an
  unknown solver type — the built-in list holds only the `dynamicMesh`
  solvers.
- ESI OpenFOAM wants `constant/turbulenceProperties`. The
  `constant/momentumTransport` spelling is OpenFOAM Foundation's and
  fails here with a missing-file error.
- Coded boundary conditions will not compile as root — see
  [`../ENVIRONMENT.md`](../ENVIRONMENT.md).
