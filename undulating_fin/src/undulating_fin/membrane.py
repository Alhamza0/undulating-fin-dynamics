"""
Elastic membrane deflection model.

Computes the deformation of the membrane connecting adjacent fin rays,
used in the elastic torque calculation.

References
----------
- Thesis Eq. (8)–(9): Elastic torque τ_{Γ,i}
- Thesis Eq. (10): Deflection function Λ(x)
"""

from __future__ import annotations

import numpy as np
from undulating_fin.params import FinParams


def membrane_deflection(D: float, x: float, q1: float, q2: float,
                        l: float, v: float) -> float:
    """Compute elastic membrane deflection Λ(x).

    Thesis Eq. (10):
        Λ(x) = (x³α)/(3D) · (1 + ν·(8 − 3α√(k²+1))/4)
              − (D²)/(2p³) · (√(k⁴+k²) − sinh⁻¹(k)) · (1 + ν·(4+α²)/4)

    where α = D/ℓ, p = √(2(1−cos(Δq))), k = (x/D)·p

    Parameters
    ----------
    D : float
        Inter-ray distance [m].
    x : float
        Evaluation height along membrane [m] (hmin or hmax).
    q1, q2 : float
        Joint angles of adjacent rays [rad].
    l : float
        Membrane natural length [m].
    v : float
        Poisson's ratio [-].

    Returns
    -------
    Lambda : float
        Membrane deflection value [m³].
    """
    alpha = D / l
    dq = q1 - q2

    p_sq = 2 * (1 - np.cos(dq))
    if p_sq < 1e-30:
        return 0.0

    p = np.sqrt(p_sq)
    kappa = (x / D) * p

    term1 = (x**3 * alpha / (3 * D)) * (
        1 + v * (8 - 3 * alpha * np.sqrt(kappa**2 + 1)) / 4
    )

    kappa_sq = kappa**2
    term2 = (D**2 / (2 * p**3)) * (
        np.sqrt(kappa_sq**2 + kappa_sq) - np.arcsinh(kappa)
    ) * (1 + v * (4 + alpha**2) / 4)

    return term1 - term2


def elastic_torques(params: FinParams, q: np.ndarray) -> np.ndarray:
    """Compute elastic membrane torques Γ(Δq) for all joints.

    Thesis Eq. (8)–(9):
        τ_{Γ,i} = Σ contributions from adjacent membranes

    Each membrane Ω_{j|j+1} exerts a torque on both joint j and joint j+1.
    The torque contribution from membrane Ω_{j|j+1} is:
        τ = (-1)^(j+1) · sin(q_j - q_{j+1}) · EK · (Λ(hmax) - Λ(hmin))
    on joint j (with sign), and the reaction on joint j+1.

    Interior rays receive contributions from TWO adjacent membranes;
    edge rays receive from only ONE.

    Parameters
    ----------
    params : FinParams
        Fin parameters.
    q : 1-D array, shape (N,)
        Joint angles [rad].

    Returns
    -------
    tau_gamma : ndarray, shape (N,)
        Elastic torques on each joint [N·m].
    """
    N = len(q)
    tau_gamma = np.zeros(N)

    # Each membrane panel Ω_{j|j+1} contributes to joints j and j+1
    for j in range(N - 1):
        dq_jk = q[j] - q[j + 1]
        delta_hmax = membrane_deflection(
            params.D, params.hmax, q[j], q[j + 1], params.l, params.v
        )
        delta_hmin = membrane_deflection(
            params.D, params.hmin, q[j], q[j + 1], params.l, params.v
        )
        membrane_torque = np.sin(dq_jk) * params.EK * (delta_hmax - delta_hmin)

        # Torque on joint j from this membrane
        tau_gamma[j] += (-1) ** (j + 1) * membrane_torque
        # Reaction torque on joint j+1 (opposite sign contribution)
        tau_gamma[j + 1] += (-1) ** (j + 2) * membrane_torque

    return tau_gamma
