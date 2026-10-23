"""
Parameter definitions for undulating fin simulations.

All parameters are organized into dataclasses with units documented.
This replaces the MATLAB global variables approach for better traceability.

References
----------
- Thesis Table 2-1: Simulation parameters
- Thesis Table 3-1: CFD simulation parameters
- Thesis Sec. 2.2.1: Fin model description
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass
class FinParams:
    """Physical and wave parameters for an undulating fin.

    Attributes
    ----------
    N : int
        Number of fin rays (links). Thesis Sec. 2.2.1.
    A : float
        Oscillation amplitude [rad]. Thesis Eq. (1).
    f : float
        Oscillation frequency [Hz]. Thesis Eq. (1).
    Phi : float
        Inter-ray phase shift [rad]. Thesis Eq. (1).
    X : float
        Angular offset (bias angle) [rad]. Thesis Eq. (1).
    D : float
        Distance between revolute joints [m]. Thesis Fig. 2-2.
    hmax : float
        Maximum membrane mounting distance from joint [m].
    hmin : float
        Minimum membrane mounting distance from joint [m].
    Cn : float
        Normal drag coefficient [-]. Thesis Sec. 3.3.4 (≈2.65).
    rho : float
        Fluid density [kg/m³].
    """

    # --- Wave parameters (Eq. 1) ---
    N: int = 3
    A: float = field(default_factory=lambda: np.deg2rad(30))
    f: float = 2.0
    Phi: float = field(default_factory=lambda: np.deg2rad(40))
    X: float = 0.0

    # --- Fin geometry (Fig. 2-2) ---
    D: float = 0.10
    hmax: float = 0.175
    hmin: float = 0.0

    # --- Hydrodynamic coefficients ---
    Cn: float = 2.65
    Cm: float = 1.0        # Added mass coefficient for Morison's eq [-]
    rho: float = 1022.0

    # --- Membrane material ---
    l: float = 0.08        # Characteristic length for DeltaFunc [m]
    v: float = 0.5         # Poisson's ratio [-]
    EK: float = 360.0      # Membrane stiffness [N/m]
    alpha: float = 1.0     # Tension level D/ℓ [-]

    # --- Fin ray mechanical properties (Sec. 2.3.1) ---
    mi: float = 0.01       # Ray mass [kg]
    li: float = 0.12       # Ray length [m]
    Izi: float = 72e-6     # Moment of inertia [kg·m²]
    bi: float = 0.02       # Viscous friction coefficient [N·m·s]
    kci: float = 0.001     # Coulomb friction coefficient [N·m]
    g: float = 9.81        # Gravitational acceleration [m/s²]

    @property
    def wavelength(self) -> float:
        """Fin wavelength λ [m]. From Eq. (3): wn = |ϕ|·(N-1)/(2π) = L/λ."""
        if abs(self.Phi) < 1e-12:
            return float("inf")
        return 2 * np.pi * self.D / abs(self.Phi)

    @property
    def wave_count(self) -> float:
        """Number of fin waves wn. Thesis Eq. (3)."""
        return abs(self.Phi) * (self.N - 1) / (2 * np.pi)

    @property
    def Mi(self) -> float:
        """Generalized inertia for one ray: mi·(li/2)² + Izi. Thesis Eq. (6)."""
        return self.mi * (self.li / 2) ** 2 + self.Izi

    @property
    def gi(self) -> float:
        """Gravity coefficient: (li/2)·mi·g. Thesis below Eq. (7)."""
        return (self.li / 2) * self.mi * self.g


@dataclass
class SimulationConfig:
    """Time-stepping configuration.

    Attributes
    ----------
    t_max : float
        Maximum simulation time [s].
    dt : float
        Time step [s].
    """

    t_max: float = 2.0
    dt: float = 0.01

    @property
    def time_vector(self) -> np.ndarray:
        """Create time vector from 0 to t_max with step dt."""
        return np.arange(0, self.t_max + self.dt / 2, self.dt)

    @property
    def n_steps(self) -> int:
        return len(self.time_vector)


@dataclass
class RobotParams:
    """Parameters for the dual-fin robot body. Thesis Ch. 6, Sec. 6.2.

    Attributes
    ----------
    R : float
        Body radius (cylinder/sphere approximation) [m].
    L : float
        Body length [m].
    m : float
        Body mass [kg].
    Cd_sphere : float
        Drag coefficient for sphere (frontal) [-].
    Cd_cylinder : float
        Drag coefficient for cylinder (lateral) [-].
    """

    R: float = 0.1
    L: float = 0.22
    m: float = 2.5
    Cd_sphere: float = 1.1
    Cd_cylinder: float = 0.70

    rho: float = 1022.0  # Fluid density, same as fin

    @property
    def added_mass_x(self) -> float:
        """Added mass in surge: (2/3)·ρ·π·R³. Thesis Sec. 6.3.1."""
        return (2 / 3) * self.rho * np.pi * self.R**3

    @property
    def added_mass_y(self) -> float:
        """Added mass in sway: ρ·π·R²·L."""
        return self.rho * np.pi * self.R**2 * self.L

    @property
    def total_mass_x(self) -> float:
        return self.m + self.added_mass_x

    @property
    def total_mass_y(self) -> float:
        return self.m + self.added_mass_y

    @property
    def total_mass_z(self) -> float:
        return self.m + self.added_mass_y  # Same as y
