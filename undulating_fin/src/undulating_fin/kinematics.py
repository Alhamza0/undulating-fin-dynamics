"""
Traveling wave kinematics for undulating fin rays.

Implements the sinusoidal motion profile for N fin rays.

References
----------
- Thesis Eq. (1): q_i(t) = A·sin(2πft − (i−1)ϕ) + X
- Thesis Eq. (2): q̇_i(t) = 2πfA·cos(2πft − (i−1)ϕ)
- Thesis Eq. (3): wn = |ϕ|·(N−1)/(2π)
"""

from __future__ import annotations

import numpy as np
from undulating_fin.params import FinParams


def compute_kinematics(
    params: FinParams,
    t: float | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute angular position, velocity, and acceleration for all rays.

    Thesis Eq. (1)–(2), extended with acceleration.

    Parameters
    ----------
    params : FinParams
        Fin parameters (N, A, f, Phi, X).
    t : float or array
        Time instant(s) [s]. Scalar or 1-D array of length nT.

    Returns
    -------
    q : ndarray, shape (nT, N)
        Joint angles [rad].
    dq : ndarray, shape (nT, N)
        Joint angular velocities [rad/s].
    ddq : ndarray, shape (nT, N)
        Joint angular accelerations [rad/s²].
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    nT = len(t)

    N = params.N
    A = params.A
    f = params.f
    Phi = params.Phi
    X = params.X

    q = np.zeros((nT, N))
    dq = np.zeros((nT, N))
    ddq = np.zeros((nT, N))

    omega = 2 * np.pi * f  # angular frequency

    for i in range(N):
        phase = i * Phi
        # Eq. (1): q_i(t) = A·sin(ωt − phase) + X
        q[:, i] = A * np.sin(omega * t - phase) + X
        # Eq. (2): q̇_i(t) = Aω·cos(ωt − phase)
        dq[:, i] = A * omega * np.cos(omega * t - phase)
        # Acceleration: q̈_i(t) = −Aω²·sin(ωt − phase)
        ddq[:, i] = -A * omega**2 * np.sin(omega * t - phase)

    return q, dq, ddq


def compute_angle_separation(q: np.ndarray) -> np.ndarray:
    """Compute angular separation Δq between adjacent rays.

    Parameters
    ----------
    q : ndarray, shape (nT, N)
        Joint angles.

    Returns
    -------
    delta_q : ndarray, shape (nT, N-1)
        Angle differences q_i − q_{i+1}.
    """
    return q[:, :-1] - q[:, 1:]
