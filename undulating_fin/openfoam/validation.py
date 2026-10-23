"""
Validation comparison: OpenFOAM results vs. analytical models.

Computes error metrics (RMSE, NRMSE, correlation) and generates
comparison plots — the same methodology used for COMSOL validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from undulating_fin.params import FinParams
from undulating_fin.kinematics import compute_kinematics
from undulating_fin.hydrodynamics import compute_total_forces, HYDRO_METHODS
from .result_parser import ForceData, load_forces, resample_to_analytical


# Force component labels matching existing CFD comparison convention
FORCE_LABELS = ("Heave", "Sway", "Surge")

# Map: analytical index → OpenFOAM index
#   analytical: [0] = Heave (sin/Fx), [1] = Sway (cos/Fy), [2] = Surge (sinΔq/Fz)
#   OpenFOAM:   [0] = Fx,             [1] = Fy,            [2] = Fz
# Default identity mapping — override via OFValidation.component_map if needed.
_DEFAULT_MAP = {0: 0, 1: 1, 2: 2}


@dataclass
class MetricSet:
    """Error metrics for one force component."""
    label: str
    rmse: float
    nrmse_pct: float    # NRMSE as percentage
    correlation: float
    max_abs_error: float


@dataclass
class ValidationResult:
    """Complete comparison result for one (case, method) pair."""
    case_name: str
    method: str
    phi_deg: float
    time: np.ndarray
    F_of: np.ndarray        # (nT, 3) OpenFOAM forces
    F_an: np.ndarray        # (nT, 3) analytical forces
    metrics: list[MetricSet]
    mean_nrmse: float


def compute_metrics(
    ref: np.ndarray,
    sig: np.ndarray,
    label: str,
) -> MetricSet:
    """Compute RMSE, NRMSE, correlation, and max error.

    Parameters
    ----------
    ref : 1-D array
        Reference signal (OpenFOAM).
    sig : 1-D array
        Signal to evaluate (analytical model).
    label : str
        Force component name.
    """
    n = min(len(ref), len(sig))
    ref, sig = ref[:n], sig[:n]

    err = sig - ref
    rmse = float(np.sqrt(np.mean(err ** 2)))
    amp = max(np.max(np.abs(ref)), 1e-12)
    nrmse = rmse / amp * 100.0

    # Pearson correlation
    if np.std(ref) < 1e-12 or np.std(sig) < 1e-12:
        corr = 0.0
    else:
        corr = float(np.corrcoef(ref, sig)[0, 1])

    return MetricSet(
        label=label,
        rmse=rmse,
        nrmse_pct=nrmse,
        correlation=corr,
        max_abs_error=float(np.max(np.abs(err))),
    )


def run_analytical(
    params: FinParams,
    time_vec: np.ndarray,
    method: str = "rft",
) -> np.ndarray:
    """Run the analytical model over *time_vec* and return (nT, 3) forces."""
    q, dq, ddq = compute_kinematics(params, time_vec)
    F = np.zeros((len(time_vec), 3))
    for i in range(len(time_vec)):
        F[i, :] = compute_total_forces(
            params, q[i, :], dq[i, :], ddq=ddq[i, :], method=method,
        )
    return F


def validate_case(
    case_dir: str | Path,
    params: FinParams,
    *,
    methods: Sequence[str] = ("rft", "morison", "lighthill"),
    component_map: dict[int, int] | None = None,
    last_n_cycles: int = 1,
) -> list[ValidationResult]:
    """Run full validation for a single OpenFOAM case.

    Parameters
    ----------
    case_dir : path-like
        Path to completed OpenFOAM case.
    params : FinParams
        Fin parameters used to generate the case.
    methods : sequence of str
        Analytical methods to compare against.
    component_map : dict, optional
        Map from analytical component index to OpenFOAM component index.
    last_n_cycles : int
        Number of trailing cycles to use for metric computation
        (avoids transient startup).

    Returns
    -------
    List of ValidationResult, one per method.
    """
    cmap = component_map or _DEFAULT_MAP
    of_data = load_forces(case_dir)

    # Build analytical time vector matching OpenFOAM time range
    dt = np.median(np.diff(of_data.time))
    time_vec = np.arange(of_data.time[0], of_data.time[-1] + dt / 2, dt)

    # Resample OpenFOAM
    F_of = resample_to_analytical(of_data, time_vec)

    # Apply component mapping
    F_of_mapped = np.column_stack([F_of[:, cmap[k]] for k in range(3)])

    # Determine last-cycle window
    spc = int(round(1.0 / (params.f * dt)))
    n = len(time_vec)
    s = max(0, n - last_n_cycles * spc)

    results = []
    case_name = Path(case_dir).name

    for method in methods:
        if method not in HYDRO_METHODS:
            continue
        F_an = run_analytical(params, time_vec, method=method)

        metrics = []
        for k, label in enumerate(FORCE_LABELS):
            m = compute_metrics(F_of_mapped[s:n, k], F_an[s:n, k], label)
            metrics.append(m)

        mean_nrmse = float(np.mean([m.nrmse_pct for m in metrics]))

        results.append(ValidationResult(
            case_name=case_name,
            method=method,
            phi_deg=float(np.rad2deg(params.Phi)),
            time=time_vec,
            F_of=F_of_mapped,
            F_an=F_an,
            metrics=metrics,
            mean_nrmse=mean_nrmse,
        ))

    return results


def batch_validate(
    cases: dict[str, tuple[Path, FinParams]],
    **kwargs,
) -> dict[str, list[ValidationResult]]:
    """Validate multiple cases.

    Parameters
    ----------
    cases : dict
        Mapping of case_name → (case_dir, FinParams).
    **kwargs
        Passed to ``validate_case``.

    Returns
    -------
    dict mapping case_name → list[ValidationResult].
    """
    all_results = {}
    for name, (case_dir, params) in cases.items():
        all_results[name] = validate_case(case_dir, params, **kwargs)
    return all_results


def results_to_dataframe(results: list[ValidationResult]):
    """Convert validation results to a pandas DataFrame for display."""
    import pandas as pd

    rows = []
    for r in results:
        for m in r.metrics:
            rows.append({
                "Case": r.case_name,
                "Φ [°]": r.phi_deg,
                "Method": r.method.upper(),
                "Component": m.label,
                "RMSE [N]": round(m.rmse, 4),
                "NRMSE [%]": round(m.nrmse_pct, 1),
                "Correlation": round(m.correlation, 4),
                "Max |Error| [N]": round(m.max_abs_error, 4),
            })
    return pd.DataFrame(rows)
