"""
Undulating Fin Hydrodynamics Package
=====================================

Kinematic and Dynamic Modeling of Biomimetic Undulating Fins
for Underwater Robots.

Based on:  Ali, A.I. (2021) "Kinematic and Dynamic Modeling of Biomimetic
           Undulating Fins for Underwater Robots", M.Sc. Thesis,
           Bahcesehir University, Dept. of Mechatronics Engineering.

Modules
-------
params       : Fin and simulation parameter definitions
kinematics   : Traveling wave formulation (Eq. 1-3)
hydrodynamics: Force and torque integration (Eq. 12-16)
dynamics     : Fin system dynamics (Eq. 4-11)
membrane     : Elastic membrane deflection (Eq. 8-10)
robot        : Dual-fin robot body dynamics (Ch. 6)
parametric   : Parametric sweep studies (Ch. 5)
plotting     : Visualization utilities
"""

from undulating_fin.params import FinParams, SimulationConfig
from undulating_fin.kinematics import compute_kinematics
from undulating_fin.hydrodynamics import (
    compute_panel_forces,
    compute_total_forces,
    compute_panel_torques,
)

__version__ = "0.1.0"

__all__ = [
    "FinParams",
    "SimulationConfig",
    "compute_kinematics",
    "compute_panel_forces",
    "compute_total_forces",
    "compute_panel_torques",
]
