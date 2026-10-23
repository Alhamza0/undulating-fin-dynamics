"""The generated case must be complete, self-consistent, and its C++
must compute the same motion as the Python reference.

A transcription that has drifted from its reference is the failure mode
that produced the previous generator's dead code, so the C++ is
extracted, compiled and run rather than eyeballed.
"""

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from undulating_fin.params import FinParams
from openfoam.case import CaseConfig, generate_case
from openfoam.geometry import membrane_displacement


PARAMS = FinParams(N=3, D=0.09, hmin=0.034, hmax=0.174,
                   A=np.deg2rad(30), f=1.0, Phi=np.deg2rad(40))


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    return generate_case(PARAMS, tmp_path_factory.mktemp("case") / "fin")


def test_all_required_files_written(case):
    expected = [
        "constant/triSurface/fin.stl",
        "constant/dynamicMeshDict",
        "constant/transportProperties",
        "constant/turbulenceProperties",
        "system/blockMeshDict",
        "system/snappyHexMeshDict",
        "system/surfaceFeatureExtractDict",
        "system/controlDict",
        "system/fvSchemes",
        "system/fvSolution",
        "system/decomposeParDict",
        "0.orig/U",
        "0.orig/p",
        "0.orig/pointDisplacement",
        "Allrun",
        "Allclean",
    ]
    missing = [f for f in expected if not (case / f).exists()]
    assert not missing, f"missing: {missing}"


def test_no_dead_solver_or_stub_motion(case):
    """Guard against the failure modes of the previous generator."""
    text = " ".join(
        (case / f).read_text()
        for f in ("system/controlDict", "0.orig/pointDisplacement",
                  "constant/dynamicMeshDict")
    )
    assert "pimpleDyMFoam" not in text          # removed from OpenFOAM
    assert "vector::zero" not in text            # the old no-op stub
    assert "placeholder" not in text.lower()
    assert "motionSolverLibs" in text            # else unknown solver type


def test_fin_patch_is_consistent_across_dictionaries(case):
    """Every reference to the fin patch must agree, or the solver dies
    at start-up on a patch that does not exist."""
    assert "patches         (fin);" in (case / "system/controlDict").read_text()
    assert "(fin)" in (case / "constant/dynamicMeshDict").read_text()
    for f in ("0.orig/U", "0.orig/p", "0.orig/pointDisplacement"):
        assert re.search(r"^\s+fin\s*$", (case / f).read_text(), re.M), f


def test_solver_entries_a_moving_mesh_needs(case):
    fv = (case / "system/fvSolution").read_text()
    assert "cellDisplacement" in fv
    assert "pcorr" in fv
    # cellDisplacement is symmetric: DILU would be rejected
    block = fv[fv.index('"cellDisplacement'):]
    assert "DIC" in block[:200] and "DILU" not in block[:200]


def test_mesh_can_resolve_the_membrane():
    """A cell coarser than the membrane cannot capture it; the config
    must refuse rather than silently produce a fin full of holes."""
    with pytest.raises(ValueError, match="cannot resolve"):
        CaseConfig(background_cell=0.04, surface_refinement=1,
                   membrane_thickness=1e-3)


def test_domain_contains_the_swept_fin(case):
    """The fin must stay inside the domain at every phase of the sweep."""
    import numpy as np

    txt = (case / "system/blockMeshDict").read_text()
    verts = [tuple(map(float, m)) for m in
             re.findall(r"\(([-\d.]+) ([-\d.]+) ([-\d.]+)\)", txt)][:8]
    lo = np.min(verts, axis=0)
    hi = np.max(verts, axis=0)

    span = (PARAMS.N - 1) * PARAMS.D
    assert lo[0] <= -PARAMS.hmax and hi[0] >= PARAMS.hmax
    assert lo[1] <= -PARAMS.hmax and hi[1] >= PARAMS.hmax
    assert lo[2] <= -span and hi[2] >= 0.0


# ── The transcription check ─────────────────────────────────────────

WSL = shutil.which("wsl")

HARNESS = r"""
#include <cmath>
#include <cstdio>
#include <cstdlib>

%s

int main()
{
    // One line of output per line of input: "t px py pz" -> "dx dy"
    double t, px, py, pz, dx, dy;
    while (scanf("%%lf %%lf %%lf %%lf", &t, &px, &py, &pz) == 4)
    {
        finMembraneDisp(t, px, py, pz, dx, dy);
        printf("%%.17g %%.17g\n", dx, dy);
    }
    return 0;
}
"""

KERNEL_BEGIN = "// --- fin motion kernel begin ---"
KERNEL_END = "// --- fin motion kernel end ---"


def _extract_kernel(case: Path) -> str:
    """Pull the C++ kernel out of the generated dictionary, undoing the
    indentation it is embedded with."""
    text = (case / "0.orig/pointDisplacement").read_text()
    body = text[text.index(KERNEL_BEGIN):text.index(KERNEL_END)]
    lines = body.splitlines()
    pad = min(
        (len(ln) - len(ln.lstrip()) for ln in lines[1:] if ln.strip()),
        default=0,
    )
    return "\n".join(
        ln[pad:] if ln.startswith(" " * pad) else ln for ln in lines
    )


def test_kernel_is_extractable_and_self_contained(case):
    """The kernel must be plain C++, with no OpenFOAM types leaking in.

    Self-containment is what allows it to be compiled and checked on
    its own, so it is worth asserting rather than assuming.
    """
    k = _extract_kernel(case)
    assert "finMembraneDisp" in k
    for foam in ("vectorField", "pointField", "labelList", "operator==",
                 "this->", "forAll", "Foam::"):
        assert foam not in k, f"kernel leaks OpenFOAM symbol: {foam!r}"


@pytest.mark.skipif(WSL is None, reason="WSL not available")
def test_cpp_transcription_matches_python(case, tmp_path):
    """Compile the generated kernel and compare it with the Python
    reference across the fin and the oscillation cycle.

    The C++ is a hand transcription of ``membrane_displacement``, and
    nothing but this test stops the two drifting apart.
    """
    cpp = tmp_path / "probe.cpp"
    cpp.write_text(HARNESS % _extract_kernel(case), encoding="utf-8", newline="\n")

    def wsl_path(path: Path) -> str:
        out = subprocess.run(
            [WSL, "-d", "Ubuntu-24.04", "--exec", "wslpath", "-a",
             str(path).replace("\\", "/")],
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()

    wp = wsl_path(cpp)
    exe = wp.rsplit("/", 1)[0] + "/probe"
    build = subprocess.run(
        [WSL, "-d", "Ubuntu-24.04", "--exec", "bash", "-lc",
         f"g++ -O2 -Wall -Wextra -o {exe} {wp} 2>&1"],
        capture_output=True, text=True,
    )
    assert build.returncode == 0, f"compile failed:\n{build.stdout}{build.stderr}"

    # Sample the membrane and the thickened solid around it, over more
    # than a full period.
    rng = np.random.default_rng(11)
    span = (PARAMS.N - 1) * PARAMS.D
    pts, times = [], []
    for t in (0.0, 0.17, 0.4, 0.63, 0.95, 1.31):
        for _ in range(40):
            z = -rng.uniform(0.0, span)
            th = rng.uniform(-np.pi, np.pi)
            r = rng.uniform(PARAMS.hmin, PARAMS.hmax)
            pts.append([r * np.cos(th), r * np.sin(th), z])
            times.append(t)
    pts = np.asarray(pts)
    times = np.asarray(times)

    # Feed the samples through a file: wsl.exe does not forward stdin
    # reliably, and a silent empty read would look like a pass.
    samples = tmp_path / "samples.txt"
    samples.write_text(
        # %.17g, not repr(): numpy 2 renders a float64 as
        # "np.float64(0.4)", which scanf cannot read — and an unparsed
        # file would come back empty rather than wrong.
        "".join(
            "%.17g %.17g %.17g %.17g\n" % (t, q[0], q[1], q[2])
            for t, q in zip(times, pts)
        ),
        encoding="utf-8", newline="\n",
    )
    run = subprocess.run(
        [WSL, "-d", "Ubuntu-24.04", "--exec", "bash", "-lc",
         f"{exe} < {wsl_path(samples)}"],
        capture_output=True, text=True, check=True,
    )

    cxx = np.array([
        [float(v) for v in line.split()]
        for line in run.stdout.strip().splitlines()
    ])
    assert len(cxx) == len(pts), "kernel did not answer every sample"

    py = np.array([
        membrane_displacement(q, PARAMS, t)[:2] for q, t in zip(pts, times)
    ])

    assert cxx == pytest.approx(py, abs=1e-13)
