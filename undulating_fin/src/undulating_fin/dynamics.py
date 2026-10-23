"""
Dynamic model of the N-ray fin system.

Implements the full dynamic equation of motion including gravity, friction,
elastic membrane effects, and hydrodynamic loads.

References
----------
- Thesis Eq. (4)–(11): Complete fin dynamics
- Thesis Eq. (7): Mq̈ + Bq̇ + G(q) + Kc·sgn(q̇) + Γ(Δq) + Φ(q,q̇) = τ
"""

from __future__ import annotations

import numpy as np
from undulating_fin.params import FinParams
from undulating_fin.membrane import elastic_torques
from undulating_fin.hydrodynamics import compute_panel_torques


def compute_mass_matrix(params: FinParams) -> np.ndarray:
    """Mass matrix M = diag(M1, M2, ..., MN). Thesis below Eq. (7).

    M_i = m_i·(l_i/2)² + I_{z,i}
    """
    return np.eye(params.N) * params.Mi


def compute_gravity_vector(params: FinParams, q: np.ndarray) -> np.ndarray:
    """Gravity vector G(q) = [g1·sin(q1), g2·sin(q2), ...]. Thesis Eq. (7)."""
    return params.gi * np.sin(q)


def compute_inertia_friction_torques(
    params: FinParams,
    q: np.ndarray,
    dq: np.ndarray,
    ddq: np.ndarray,
) -> np.ndarray:
    """Compute L(q̈,q̇,q) = Mq̈ + Bq̇ + Kc·sgn(q̇) + G(q).

    Thesis Eq. (11): The force on each joint due to gravity, inertia, friction.

    Parameters
    ----------
    params : FinParams
    q, dq, ddq : 1-D arrays, shape (N,)
        Joint angles, velocities, accelerations.

    Returns
    -------
    L : ndarray, shape (N,)
        Torques due to inertia, friction, and gravity [N·m].
    """
    N = len(q)
    M = params.Mi * ddq                        # M·q̈
    B = params.bi * dq                          # B·q̇
    Kc = params.kci * np.sign(dq)               # Kc·sgn(q̇)
    G = compute_gravity_vector(params, q)

    return M + B + Kc + G


def compute_total_torques(
    params: FinParams,
    q: np.ndarray,
    dq: np.ndarray,
    ddq: np.ndarray,
    include_hydrodynamic: bool = True,
) -> np.ndarray:
    """Compute total joint torques τ = L + Γ + Φ.

    Thesis Eq. (7)/(11):
        τ = L(q̈,q̇,q) + Γ(Δq) + Φ(q,q̇)

    Parameters
    ----------
    params : FinParams
    q, dq, ddq : 1-D arrays, shape (N,)
    include_hydrodynamic : bool
        If True, compute hydrodynamic torques (expensive). Default True.

    Returns
    -------
    tau : ndarray, shape (N,)
        Total required actuation torques [N·m].
    """
    N = len(q)

    # L: inertia + friction + gravity
    L = compute_inertia_friction_torques(params, q, dq, ddq)

    # Γ: elastic membrane torques
    Gamma = elastic_torques(params, q)

    # Φ: hydrodynamic torques
    Phi = np.zeros(N)
    if include_hydrodynamic and N >= 2:
        for i in range(N - 1):
            tau1, tau2 = compute_panel_torques(
                params.D, params.hmin, params.hmax,
                q[i], q[i + 1], dq[i], dq[i + 1],
                params.Cn, params.rho,
            )
            Phi[i] += tau1
            Phi[i + 1] += tau2

    return L + Gamma + Phi


def run_time_domain(
    params: FinParams,
    time_vec: np.ndarray,
    *,
    method: str = "rft",
) -> dict[str, np.ndarray]:
    """Run a complete time-domain simulation.

    Computes kinematics, forces, and torques over time.

    Parameters
    ----------
    params : FinParams
    time_vec : 1-D array
        Time values [s].
    method : str
        Hydrodynamic model: 'rft', 'morison', or 'lighthill'.

    Returns
    -------
    results : dict
        Keys: 'time', 'q', 'dq', 'ddq', 'F_N', 'T_phi', 'T_gamma', 'L', 'tau', 'method'
    """
    from undulating_fin.kinematics import compute_kinematics
    from undulating_fin.hydrodynamics import compute_total_forces

    nT = len(time_vec)
    N = params.N

    # Compute full kinematics
    q, dq, ddq = compute_kinematics(params, time_vec)

    # Allocate results
    F_N = np.zeros((nT, 3))
    T_phi = np.zeros((nT, N))
    T_gamma = np.zeros((nT, N))
    L_vec = np.zeros((nT, N))
    tau_vec = np.zeros((nT, N))

    for t_idx in range(nT):
        qt = q[t_idx, :]
        dqt = dq[t_idx, :]
        ddqt = ddq[t_idx, :]

        # Hydrodynamic forces
        F_N[t_idx, :] = compute_total_forces(
            params, qt, dqt, ddq=ddqt, method=method,
        )

        # Hydrodynamic torques
        for i in range(N - 1):
            tau1, tau2 = compute_panel_torques(
                params.D, params.hmin, params.hmax,
                qt[i], qt[i + 1], dqt[i], dqt[i + 1],
                params.Cn, params.rho,
            )
            T_phi[t_idx, i] += tau1
            T_phi[t_idx, min(i + 1, N - 1)] += tau2

        # Elastic torques
        T_gamma[t_idx, :] = elastic_torques(params, qt)

        # Inertia/friction/gravity
        L_vec[t_idx, :] = compute_inertia_friction_torques(params, qt, dqt, ddqt)

        # Total torque
        tau_vec[t_idx, :] = L_vec[t_idx, :] + T_gamma[t_idx, :] + T_phi[t_idx, :]

    return {
        "time": time_vec,
        "q": q,
        "dq": dq,
        "ddq": ddq,
        "F_N": F_N,
        "T_phi": T_phi,
        "T_gamma": T_gamma,
        "L": L_vec,
        "tau": tau_vec,
        "method": method,
    }
