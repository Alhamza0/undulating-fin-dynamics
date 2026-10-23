"""
Parametric sweep studies for undulating fin optimization.

Provides functions to sweep over fin parameters (Φ, A, f, X, N, D, H)
and compute maximum hydrodynamic forces — replacing the 8+ separate
MATLAB functions (PhiFN.m, AmpFN.m, freqFN.m, etc.).

References
----------
- Thesis Ch. 5: Parametric sweep
- Thesis Table 5-1: Effect of parameters on forces
- Thesis Sec. 5.1–5.3: Undulating parameters, membrane, ray count
"""

from __future__ import annotations

import copy
from dataclasses import replace
from typing import Callable

import numpy as np
from undulating_fin.params import FinParams, SimulationConfig
from undulating_fin.kinematics import compute_kinematics
from undulating_fin.hydrodynamics import compute_total_forces


def sweep_1d(
    base_params: FinParams,
    sim_config: SimulationConfig,
    param_name: str,
    values: np.ndarray,
    verbose: bool = True,
    *,
    method: str = "rft",
) -> dict[str, np.ndarray]:
    """Perform a 1-D parametric sweep over a single parameter.

    For each value in *values*, sets ``base_params.<param_name>`` to that
    value, runs a time-domain force computation, and records the maximum
    of each force component.

    Parameters
    ----------
    base_params : FinParams
        Base parameters (all other params stay fixed).
    sim_config : SimulationConfig
        Time settings.
    param_name : str
        Name of the FinParams attribute to sweep (e.g. 'Phi', 'A', 'f').
    values : 1-D array
        Values to sweep over.
    verbose : bool
        Print progress.

    Returns
    -------
    results : dict
        'values' : sweep values
        'F_max'  : (n, 3) max [Heave, Sway, Surge] at each value
        'param_name' : swept parameter name

    Examples
    --------
    >>> # Equivalent to MATLAB PhiFN.m:
    >>> r = sweep_1d(params, sim, 'Phi', np.linspace(0, np.deg2rad(80), 100))
    >>> plt.plot(np.rad2deg(r['values']), r['F_max'][:, 2])  # Surge vs Phi
    """
    time_vec = sim_config.time_vector
    nT = sim_config.n_steps
    n_vals = len(values)
    F_max = np.zeros((n_vals, 3))

    for idx, val in enumerate(values):
        if verbose and idx % max(1, n_vals // 10) == 0:
            print(f"  sweep {param_name}: {idx+1}/{n_vals}")

        # Create modified params
        p = copy.copy(base_params)
        setattr(p, param_name, val)

        # Compute kinematics for all time steps
        q, dq, ddq = compute_kinematics(p, time_vec)

        # Compute forces at each time step
        F_time = np.zeros((nT, 3))
        for t_idx in range(nT):
            F_time[t_idx, :] = compute_total_forces(
                p, q[t_idx, :], dq[t_idx, :],
                ddq=ddq[t_idx, :], method=method,
            )

        F_max[idx, :] = np.max(F_time, axis=0)

    return {
        "values": values,
        "F_max": F_max,
        "param_name": param_name,
    }


def sweep_2d(
    base_params: FinParams,
    sim_config: SimulationConfig,
    param1_name: str,
    values1: np.ndarray,
    param2_name: str,
    values2: np.ndarray,
    force_component: int = 2,
    verbose: bool = True,
    *,
    method: str = "rft",
) -> dict[str, np.ndarray]:
    """Perform a 2-D parametric sweep over two parameters.

    Replaces MATLAB Phi_A_FUN.m, Phi_D_FUN.m, Phi_F_FUN.m, Phi_H_FUN.m.

    Parameters
    ----------
    base_params : FinParams
    sim_config : SimulationConfig
    param1_name, param2_name : str
        Names of parameters to sweep.
    values1, values2 : 1-D arrays
        Sweep values for each parameter.
    force_component : int
        0=Heave, 1=Sway, 2=Surge (default).
    verbose : bool

    Returns
    -------
    results : dict
        'values1', 'values2' : sweep grids
        'F_2d' : (n1, n2) max force values
        'param1_name', 'param2_name' : parameter names

    Examples
    --------
    >>> # Equivalent to Phi_A_FUN.m:
    >>> r = sweep_2d(p, s, 'Phi', phi_vec, 'A', a_vec, force_component=2)
    """
    time_vec = sim_config.time_vector
    nT = sim_config.n_steps
    n1, n2 = len(values1), len(values2)
    F_2d = np.zeros((n1, n2))

    for j, v2 in enumerate(values2):
        if verbose:
            print(f"  sweep {param2_name}[{j+1}/{n2}]")
        for i, v1 in enumerate(values1):
            p = copy.copy(base_params)
            setattr(p, param1_name, v1)
            setattr(p, param2_name, v2)

            q, dq, ddq = compute_kinematics(p, time_vec)
            F_time = np.zeros((nT, 3))
            for t_idx in range(nT):
                F_time[t_idx, :] = compute_total_forces(
                    p, q[t_idx, :], dq[t_idx, :],
                    ddq=ddq[t_idx, :], method=method,
                )

            F_2d[i, j] = np.max(F_time[:, force_component])

    return {
        "values1": values1,
        "values2": values2,
        "F_2d": F_2d,
        "param1_name": param1_name,
        "param2_name": param2_name,
    }


def sweep_ray_count(
    base_params: FinParams,
    sim_config: SimulationConfig,
    N_values: list[int] | np.ndarray,
    phi_values: np.ndarray,
    verbose: bool = True,
    *,
    method: str = "rft",
) -> dict[str, np.ndarray]:
    """Sweep over number of rays N × phase shift Φ.

    Replaces MATLAB Number_of_Rays.m. Computes max forces for each (N, Φ).

    Parameters
    ----------
    base_params : FinParams
    sim_config : SimulationConfig
    N_values : list of int
        Ray counts to test (e.g. [2, 3, 4, 5, 6]).
    phi_values : 1-D array
        Phase shift values to sweep.

    Returns
    -------
    results : dict with keys per force component per N.
    """
    time_vec = sim_config.time_vector
    nT = sim_config.n_steps
    n_phi = len(phi_values)
    n_N = len(N_values)

    # F_heave[i_phi, i_N], F_sway[...], F_surge[...]
    F_heave = np.zeros((n_phi, n_N))
    F_sway = np.zeros((n_phi, n_N))
    F_surge = np.zeros((n_phi, n_N))

    for j, N in enumerate(N_values):
        if verbose:
            print(f"  N={N} ({j+1}/{n_N})")
        for i, phi in enumerate(phi_values):
            p = copy.copy(base_params)
            p.N = N
            p.Phi = phi

            q, dq, ddq = compute_kinematics(p, time_vec)
            F_time = np.zeros((nT, 3))
            for t_idx in range(nT):
                F_time[t_idx, :] = compute_total_forces(
                    p, q[t_idx, :], dq[t_idx, :],
                    ddq=ddq[t_idx, :], method=method,
                )

            F_heave[i, j] = np.max(F_time[:, 0])
            F_sway[i, j] = np.max(F_time[:, 1])
            F_surge[i, j] = np.max(F_time[:, 2])

    return {
        "phi_values": phi_values,
        "N_values": np.array(N_values),
        "F_heave": F_heave,
        "F_sway": F_sway,
        "F_surge": F_surge,
    }


# ════════════════════════════════════════════════════════════════════
# Force stability / smoothness metrics
# ════════════════════════════════════════════════════════════════════

def _force_time_history(
    params: FinParams,
    sim_config: SimulationConfig,
    method: str = "rft",
) -> np.ndarray:
    """Return full force time-history (nT, 3) for the given params."""
    t = sim_config.time_vector
    q, dq, ddq = compute_kinematics(params, t)
    nT = len(t)
    F = np.zeros((nT, 3))
    for ti in range(nT):
        F[ti, :] = compute_total_forces(
            params, q[ti, :], dq[ti, :],
            ddq=ddq[ti, :], method=method,
        )
    return F


def compute_stability_metrics(F_time: np.ndarray, dt: float) -> dict:
    """Compute force stability / smoothness metrics from a time-history.

    Parameters
    ----------
    F_time : (nT, 3) array
        Force time-history [Heave, Sway, Surge].
    dt : float
        Time step [s].

    Returns
    -------
    dict with keys per component ('Heave', 'Sway', 'Surge'):
        'mean'       : mean force [N]
        'std'        : standard deviation [N]
        'peak2peak'  : max − min [N]
        'ripple_pct' : 100 × std / |mean| (∞ if mean ≈ 0)
        'smoothness' : RMS of dF/dt (lower = smoother)
        'sign_changes' : number of zero-crossings (force reversals)
        'pos_fraction' : fraction of time the force is positive
    """
    labels = ["Heave", "Sway", "Surge"]
    metrics = {}
    for k, lbl in enumerate(labels):
        f = F_time[:, k]
        mean = float(np.mean(f))
        std = float(np.std(f))
        p2p = float(np.max(f) - np.min(f))

        # Ripple: how much ± oscillation relative to mean
        ripple = 100.0 * std / abs(mean) if abs(mean) > 1e-12 else float("inf")

        # Smoothness: RMS of force derivative (jerk-related)
        df = np.diff(f) / dt
        smoothness = float(np.sqrt(np.mean(df ** 2)))

        # Sign changes (zero-crossings) — force direction reversals
        signs = np.sign(f)
        sign_changes = int(np.sum(np.abs(np.diff(signs)) > 0))

        # Fraction of time force is positive (unidirectional = 1.0)
        pos_frac = float(np.sum(f > 0) / len(f))

        metrics[lbl] = {
            "mean": mean,
            "std": std,
            "peak2peak": p2p,
            "ripple_pct": ripple,
            "smoothness": smoothness,
            "sign_changes": sign_changes,
            "pos_fraction": pos_frac,
        }
    return metrics


def sweep_N_stability(
    base_params: FinParams,
    sim_config: SimulationConfig,
    N_values: list[int] | np.ndarray,
    verbose: bool = True,
    *,
    method: str = "rft",
) -> dict:
    """Sweep ray count N and compute force stability metrics.

    This answers: as N ↑ (wave becomes more continuous / sine-like),
    do the forces become smoother and more unidirectional?

    Parameters
    ----------
    base_params : FinParams
    sim_config : SimulationConfig
    N_values : list of int
        Ray counts to evaluate.
    method : str
        Hydrodynamic model ('rft', 'morison', 'lighthill').

    Returns
    -------
    dict with:
        'N_values' : array of N
        'metrics'  : list of stability dicts (one per N)
        'F_series' : list of (nT, 3) force arrays (one per N)
        'wave_counts' : array — number of waves on fin for each N
    """
    dt = sim_config.dt
    results_metrics = []
    F_list = []
    wn_list = []

    for N in N_values:
        if verbose:
            print(f"  N={N} rays...")
        p = copy.copy(base_params)
        p.N = int(N)

        F = _force_time_history(p, sim_config, method=method)
        m = compute_stability_metrics(F, dt)

        results_metrics.append(m)
        F_list.append(F)
        wn_list.append(p.wave_count)

    return {
        "N_values": np.array(N_values),
        "metrics": results_metrics,
        "F_series": F_list,
        "wave_counts": np.array(wn_list),
    }


def sweep_geometry_stability(
    base_params: FinParams,
    sim_config: SimulationConfig,
    param_name: str,
    values: np.ndarray,
    verbose: bool = True,
    *,
    method: str = "rft",
) -> dict:
    """Sweep a physical dimension (D, hmax, hmin, li) and compute stability.

    Parameters
    ----------
    param_name : str
        FinParams attribute to sweep ('D', 'hmax', 'hmin', 'li').
    values : array
        Values to sweep (in SI units).

    Returns
    -------
    dict with:
        'param_name', 'values' : swept parameter info
        'metrics' : list of stability dicts
        'F_series' : list of (nT, 3) force arrays
    """
    dt = sim_config.dt
    results_metrics = []
    F_list = []

    for idx, val in enumerate(values):
        if verbose and idx % max(1, len(values) // 10) == 0:
            print(f"  {param_name}={val:.4f} ({idx+1}/{len(values)})")
        p = copy.copy(base_params)
        setattr(p, param_name, val)

        F = _force_time_history(p, sim_config, method=method)
        m = compute_stability_metrics(F, dt)

        results_metrics.append(m)
        F_list.append(F)

    return {
        "param_name": param_name,
        "values": values,
        "metrics": results_metrics,
        "F_series": F_list,
    }
