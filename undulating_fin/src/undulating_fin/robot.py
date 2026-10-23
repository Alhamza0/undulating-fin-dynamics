"""
Dual-fin robot body dynamics.

Simulates the 3-DOF translational motion of an underwater robot
propelled by two undulating fins (right and left).

References
----------
- Thesis Ch. 6: Motion Control
- Thesis Sec. 6.2: Mini-robot model
- Thesis Sec. 6.3: Basic motion and its control
- Thesis Eq. (18)–(19): Right/left fin motion
- Thesis Table 6-1: Motion parameter settings
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from undulating_fin.params import FinParams, RobotParams, SimulationConfig
from undulating_fin.kinematics import compute_kinematics
from undulating_fin.hydrodynamics import compute_total_forces


@dataclass
class MotionConfig:
    """Configuration for a specific robot maneuver.

    Thesis Table 6-1 defines parameter settings for each basic motion.
    """
    right_fin: FinParams
    left_fin: FinParams
    description: str = ""

    @staticmethod
    def forward(
        A: float = np.deg2rad(30),
        f: float = 1.0,
        Phi: float = np.deg2rad(60),
        **kwargs,
    ) -> "MotionConfig":
        """Forward motion: both fins undulate identically with ϕ > 0."""
        right = FinParams(A=A, f=f, Phi=Phi, X=0.0, **kwargs)
        left = FinParams(A=A, f=f, Phi=Phi, X=0.0, **kwargs)
        return MotionConfig(right, left, "Forward")

    @staticmethod
    def backward(
        A: float = np.deg2rad(30),
        f: float = 1.0,
        Phi: float = np.deg2rad(60),
        **kwargs,
    ) -> "MotionConfig":
        """Backward motion: both fins with ϕ < 0."""
        right = FinParams(A=A, f=f, Phi=-Phi, X=0.0, **kwargs)
        left = FinParams(A=A, f=f, Phi=-Phi, X=0.0, **kwargs)
        return MotionConfig(right, left, "Backward")

    @staticmethod
    def rotate_cw(
        A: float = np.deg2rad(30),
        f: float = 1.0,
        Phi: float = np.deg2rad(60),
        **kwargs,
    ) -> "MotionConfig":
        """Clockwise in-place rotation: ϕ_R = +Θ, ϕ_L = −Θ."""
        right = FinParams(A=A, f=f, Phi=Phi, X=0.0, **kwargs)
        left = FinParams(A=A, f=f, Phi=-Phi, X=0.0, **kwargs)
        return MotionConfig(right, left, "Clockwise rotation")


def simulate_robot(
    robot: RobotParams,
    motion: MotionConfig,
    sim: SimulationConfig,
) -> dict[str, np.ndarray]:
    """Simulate robot body motion driven by dual undulating fins.

    Uses forward Euler integration (as in original MATLAB positionandVelocity.m).
    Robot body modeled as cylinder with added mass and quadratic drag.

    Equations of motion (Thesis Sec. 6.3.1):
        m_x · ẍ = F_surge_R + F_surge_L − f_drag,x
        m_y · ÿ = F_heave_R − F_heave_L − f_drag,y
        m_z · z̈ = F_sway_R + F_sway_L

    Parameters
    ----------
    robot : RobotParams
        Robot body parameters.
    motion : MotionConfig
        Right and left fin motion settings.
    sim : SimulationConfig
        Time parameters.

    Returns
    -------
    results : dict
        Keys: 'time', 'pos' (nT,3), 'vel' (nT,3), 'acc' (nT,3),
              'F_right' (nT,3), 'F_left' (nT,3), 'F_total' (nT,3)
    """
    time_vec = sim.time_vector
    nT = sim.n_steps
    dt = sim.dt

    # Compute fin forces for right and left fins
    q_R, dq_R, _ = compute_kinematics(motion.right_fin, time_vec)
    q_L, dq_L, _ = compute_kinematics(motion.left_fin, time_vec)

    F_right = np.zeros((nT, 3))
    F_left = np.zeros((nT, 3))

    for t_idx in range(nT):
        F_right[t_idx, :] = compute_total_forces(
            motion.right_fin, q_R[t_idx, :], dq_R[t_idx, :]
        )
        F_left[t_idx, :] = compute_total_forces(
            motion.left_fin, q_L[t_idx, :], dq_L[t_idx, :]
        )

    # Robot body dynamics — forward Euler
    # Axis mapping (Thesis Fig. 6-1):
    #   Surge (x): F_right[surge] + F_left[surge]       — index 2
    #   Sway  (y): F_right[heave] - F_left[heave]       — index 0
    #   Heave (z): F_right[sway]  + F_left[sway]        — index 1
    FX = F_right[:, 2] + F_left[:, 2]   # Surge  (fin F[2])
    FY = F_right[:, 0] - F_left[:, 0]   # Sway   (differential fin heave, F[0])
    FZ = F_right[:, 1] + F_left[:, 1]   # Heave  (fin sway, F[1])

    mx = robot.total_mass_x
    my = robot.total_mass_y
    mz = robot.total_mass_z

    pos = np.zeros((nT, 3))
    vel = np.zeros((nT, 3))
    acc = np.zeros((nT, 3))

    # Initial acceleration
    acc[0, 0] = FX[0] / mx
    acc[0, 1] = FY[0] / my
    acc[0, 2] = FZ[0] / mz

    for n in range(1, nT):
        # Update position
        pos[n, :] = pos[n - 1, :] + vel[n - 1, :] * dt

        # Update velocity
        vel[n, :] = vel[n - 1, :] + acc[n - 1, :] * dt

        # Drag forces
        fd_x = 0.5 * robot.rho * vel[n - 1, 0] ** 2 * robot.Cd_sphere * np.pi * robot.R**2
        fd_y = 0.5 * robot.rho * vel[n - 1, 1] ** 2 * robot.Cd_cylinder * 2 * np.pi * robot.R * robot.L
        fd_z = 0.5 * robot.rho * vel[n - 1, 2] ** 2 * robot.Cd_cylinder * 2 * np.pi * robot.R * robot.L

        # Update acceleration
        acc[n, 0] = (FX[n] - fd_x) / mx
        acc[n, 1] = (FY[n] - fd_y) / my
        acc[n, 2] = FZ[n] / mz

    return {
        "time": time_vec,
        "pos": pos,
        "vel": vel,
        "acc": acc,
        "F_right": F_right,
        "F_left": F_left,
        "F_total": np.column_stack([FX, FY, FZ]),
    }
