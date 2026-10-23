"""
Example: Time-domain analysis of a 3-ray undulating fin.

Reproduces Thesis Sec. 2.4 (Fig. 2-5, 2-6) — hydrodynamic forces
and torques for a two-ray fin with and without phase shift.

Usage
-----
    python -m examples.run_analysis
    # or
    python examples/run_analysis.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import matplotlib.pyplot as plt
from undulating_fin.params import FinParams, SimulationConfig
from undulating_fin.dynamics import run_time_domain
from undulating_fin.plotting import plot_time_domain


def main():
    # ===== Case 1: ϕ = 0° (no fin wave, plate oscillation) =====
    # Thesis Fig. 2-5: A=30°, ϕ=0°, D=9cm, f=1Hz, X=0°, N=2
    print("=" * 60)
    print("Case 1: Two-ray fin, ϕ = 0° (no undulation)")
    print("=" * 60)

    params_phi0 = FinParams(
        N=2,
        A=np.deg2rad(30),
        f=1.0,
        Phi=np.deg2rad(0),
        X=0.0,
        D=0.09,
        hmax=0.174,
        hmin=0.034,
        Cn=2.65,
        rho=1000.0,
    )
    sim = SimulationConfig(t_max=3.0, dt=0.01)

    print("Running simulation (this may take a few minutes)...")
    results_phi0 = run_time_domain(params_phi0, sim.time_vector)

    fig1 = plot_time_domain(
        results_phi0,
        title=r"$A=30°,\ \phi=0°,\ D=9\mathrm{cm},\ f=1\mathrm{Hz}$ (Fig. 2-5)"
    )

    # ===== Case 2: ϕ = 50° (undulating fin) =====
    # Thesis Fig. 2-6: A=30°, ϕ=50°, D=9cm, f=1Hz, X=0°, N=2
    print()
    print("=" * 60)
    print("Case 2: Two-ray fin, ϕ = 50° (undulating)")
    print("=" * 60)

    params_phi50 = FinParams(
        N=2,
        A=np.deg2rad(30),
        f=1.0,
        Phi=np.deg2rad(50),
        X=0.0,
        D=0.09,
        hmax=0.174,
        hmin=0.034,
        Cn=2.65,
        rho=1000.0,
    )

    print("Running simulation...")
    results_phi50 = run_time_domain(params_phi50, sim.time_vector)

    fig2 = plot_time_domain(
        results_phi50,
        title=r"$A=30°,\ \phi=50°,\ D=9\mathrm{cm},\ f=1\mathrm{Hz}$ (Fig. 2-6)"
    )

    # ===== Print key results =====
    F = results_phi50["F_N"]
    print(f"\nMax Sway  force: {np.max(F[:, 0]):.4f} N")
    print(f"Max Heave force: {np.max(F[:, 1]):.4f} N")
    print(f"Max Surge force: {np.max(F[:, 2]):.4f} N")
    print(f"Mean Surge force: {np.mean(F[:, 2]):.4f} N")

    plt.show()


if __name__ == "__main__":
    main()
