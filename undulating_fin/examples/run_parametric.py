"""
Example: Parametric sweep studies.

Reproduces Thesis Ch. 5 — effect of undulating parameters on forces.

Usage
-----
    python examples/run_parametric.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import matplotlib.pyplot as plt
from undulating_fin.params import FinParams, SimulationConfig
from undulating_fin.parametric import sweep_1d, sweep_ray_count
from undulating_fin.plotting import plot_parametric_1d, plot_ray_count_comparison


def main():
    # Base parameters (Thesis Table 2-1)
    params = FinParams(
        N=3,
        A=np.deg2rad(30),
        f=1.0,
        Phi=np.deg2rad(40),
        X=0.0,
        D=0.10,
        hmax=0.175,
        hmin=0.0,
        Cn=2.65,
        rho=1000.0,
    )
    sim = SimulationConfig(t_max=2.0, dt=0.01)

    # ===== 1. Phase shift sweep (Fig. 5-5) =====
    print("Sweeping phase shift ϕ (0° → 80°)...")
    phi_vec = np.linspace(0, np.deg2rad(80), 30)
    r_phi = sweep_1d(params, sim, "Phi", phi_vec)
    fig1 = plot_parametric_1d(r_phi, xlabel=r"$\phi$ [°]", convert_deg=True)
    fig1.suptitle("Phase shift effect (Fig. 5-5)")

    # ===== 2. Frequency sweep (Fig. 5-3a) =====
    print("Sweeping frequency f (0.5 → 4 Hz)...")
    freq_vec = np.linspace(0.5, 4, 20)
    r_freq = sweep_1d(params, sim, "f", freq_vec)
    fig2 = plot_parametric_1d(r_freq, xlabel="$f$ [Hz]")
    fig2.suptitle("Frequency effect (Fig. 5-3a)")

    # ===== 3. Amplitude sweep (Fig. 5-3b) =====
    print("Sweeping amplitude A (0° → 40°)...")
    amp_vec = np.linspace(np.deg2rad(5), np.deg2rad(40), 20)
    r_amp = sweep_1d(params, sim, "A", amp_vec)
    fig3 = plot_parametric_1d(r_amp, xlabel="$A$ [°]", convert_deg=True)
    fig3.suptitle("Amplitude effect (Fig. 5-3b)")

    # ===== 4. Ray count comparison (Fig. 5-15) =====
    print("Sweeping ray count N × ϕ...")
    phi_vec_N = np.linspace(0, np.deg2rad(80), 20)
    r_N = sweep_ray_count(params, sim, N_values=[2, 3, 4, 5, 6], phi_values=phi_vec_N)
    fig4 = plot_ray_count_comparison(r_N)

    plt.show()


if __name__ == "__main__":
    main()
