r"""
Compare OpenFOAM forces against the analytical models.

The point of the CFD campaign is a like-for-like comparison, and the
two sides do not naturally speak the same language:

* OpenFOAM writes a force every time step from the very first, when the
  fluid is still at rest and the impulsive start is ringing through the
  domain.  Only the final cycle is a periodic state worth comparing.
* The analytical model is evaluated at whatever instants we ask for, so
  it has to be sampled on the CFD's own time base.

Both sides report the force **on the fin**: OpenFOAM's ``forces``
function object assembles ``rho*Sf*(p - pRef)`` with ``Sf`` pointing out
of the fluid and into the wall, and thesis Eq. (12) has
``f_n = -1/2 rho Cd en vn |vn|``, which is invariant under flipping
``en``.  So the comparison needs no sign convention flag; a sign
disagreement is a real disagreement.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from undulating_fin.hydrodynamics import compute_total_forces
from undulating_fin.kinematics import compute_kinematics
from undulating_fin.params import FinParams

__all__ = [
    "ForceSeries",
    "load_openfoam_forces",
    "analytical_forces",
    "last_cycle",
    "metrics",
    "compare",
]

#: Force component order used throughout the project.
COMPONENTS = ("heave", "sway", "surge")   # x, y, z


@dataclass
class ForceSeries:
    """Force history in newtons, ``F`` of shape ``(n, 3)`` as (x, y, z)."""

    t: np.ndarray
    F: np.ndarray

    def __post_init__(self) -> None:
        self.t = np.asarray(self.t, dtype=float)
        self.F = np.asarray(self.F, dtype=float)
        if self.F.shape != (self.t.size, 3):
            raise ValueError(
                f"F must be (n, 3) matching t of length {self.t.size}, "
                f"got {self.F.shape}"
            )

    def __len__(self) -> int:
        return self.t.size


# ── Loading ─────────────────────────────────────────────────────────

def load_openfoam_forces(case: str | Path) -> ForceSeries:
    """Read ``postProcessing/forces/*/force.dat`` from a case.

    Restarted runs leave several time directories; they are
    concatenated in time order and duplicated instants dropped, keeping
    the later value.
    """
    case = Path(case)
    files = sorted((case / "postProcessing/forces").glob("*/force.dat"))
    if not files:
        raise FileNotFoundError(f"no force.dat under {case}/postProcessing/forces")

    chunks = []
    for f in files:
        d = np.loadtxt(f, comments="#")
        if d.ndim == 1:
            d = d[None, :]
        chunks.append(d[:, :4])

    d = np.concatenate(chunks, axis=0)
    order = np.argsort(d[:, 0], kind="stable")
    d = d[order]
    keep = np.append(np.diff(d[:, 0]) > 0, True)   # last wins on duplicates
    d = d[keep]
    return ForceSeries(d[:, 0], d[:, 1:4])


def analytical_forces(
    params: FinParams, t: np.ndarray, method: str = "morison",
) -> ForceSeries:
    """Evaluate an analytical model on the given time base."""
    t = np.asarray(t, dtype=float)
    q, dq, ddq = compute_kinematics(params, t)
    # ddq is not optional: without it the added-mass term vanishes and
    # 'morison' silently collapses onto 'rft'.
    F = np.array([
        compute_total_forces(params, q[k], dq[k], ddq=ddq[k], method=method)
        for k in range(t.size)
    ])
    return ForceSeries(t, F)


# ── Windowing ───────────────────────────────────────────────────────

def last_cycle(series: ForceSeries, f: float, n_cycles: int = 1) -> ForceSeries:
    """Restrict to the final ``n_cycles`` oscillation periods.

    Everything before that is contaminated by the impulsive start: the
    fluid begins at rest and the first time steps of the probe runs
    carried forces two orders of magnitude above the periodic state.
    """
    if f <= 0:
        raise ValueError("f must be positive")
    t_end = series.t[-1]
    t0 = t_end - n_cycles / f
    if t0 < series.t[0] - 1e-12:
        raise ValueError(
            f"run covers {series.t[0]:.3f}–{t_end:.3f} s, too short for "
            f"{n_cycles} cycle(s) at f = {f} Hz"
        )
    m = series.t >= t0
    return ForceSeries(series.t[m], series.F[m])


def is_periodic(series: ForceSeries, f: float, rtol: float = 0.05) -> bool:
    """True if the last two cycles agree to within *rtol*.

    The honest test that the start-up transient has died away.  Without
    it, a comparison may be measuring how the fluid was set in motion
    rather than how the fin drives it.
    """
    try:
        a = last_cycle(series, f, 2)
    except ValueError:
        return False
    t_mid = a.t[-1] - 1.0 / f
    first, second = a.t < t_mid, a.t >= t_mid
    if first.sum() < 8 or second.sum() < 8:
        return False

    grid = np.linspace(0.0, 1.0, 64, endpoint=False)
    out = []
    for m in (first, second):
        ph = (a.t[m] * f) % 1.0
        o = np.argsort(ph)
        out.append(np.stack([
            np.interp(grid, ph[o], a.F[m][o, c]) for c in range(3)
        ], axis=1))

    scale = np.abs(out[0]).max()
    if scale == 0.0:
        return True
    return bool(np.abs(out[0] - out[1]).max() / scale <= rtol)


# ── Metrics ─────────────────────────────────────────────────────────

def metrics(cfd: ForceSeries, model: ForceSeries) -> dict[str, dict[str, float]]:
    """Per-component agreement, in the same terms the paper reports.

    NRMSE is normalised by the CFD peak-to-peak range, the same
    normalisation used for the published COMSOL comparison, so the
    numbers stay comparable with it.
    """
    if not np.allclose(cfd.t, model.t):
        raise ValueError("series must share a time base; sample the model on cfd.t")

    out: dict[str, dict[str, float]] = {}
    for i, name in enumerate(COMPONENTS):
        a, b = cfd.F[:, i], model.F[:, i]
        err = b - a
        rng = float(a.max() - a.min())
        rmse = float(np.sqrt(np.mean(err ** 2)))

        if np.std(a) < 1e-14 or np.std(b) < 1e-14:
            corr = float("nan")
        else:
            corr = float(np.corrcoef(a, b)[0, 1])

        out[name] = {
            "rmse_N": rmse,
            "nrmse_pct": float("nan") if rng == 0 else 100.0 * rmse / rng,
            "max_abs_error_N": float(np.abs(err).max()),
            "correlation": corr,
            "mean_cfd_N": float(a.mean()),
            "mean_model_N": float(b.mean()),
            "amp_cfd_N": rng,
            "amp_model_N": float(b.max() - b.min()),
        }
    return out


def compare(
    case: str | Path,
    params: FinParams,
    methods: tuple[str, ...] = ("rft", "morison", "lighthill"),
    n_cycles: int = 1,
) -> dict:
    """Load a case, take its final cycle, and score each model on it."""
    cfd_all = load_openfoam_forces(case)
    cfd = last_cycle(cfd_all, params.f, n_cycles)

    result = {
        "case": str(case),
        "t_range": (float(cfd_all.t[0]), float(cfd_all.t[-1])),
        "n_samples": len(cfd),
        "periodic": is_periodic(cfd_all, params.f),
        "models": {},
    }
    for m in methods:
        result["models"][m] = metrics(cfd, analytical_forces(params, cfd.t, m))
    return result
