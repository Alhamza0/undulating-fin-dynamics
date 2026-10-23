"""
Visualization utilities for undulating fin simulations.

Provides standardized plot functions for thesis-quality figures.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


def plot_time_domain(results: dict, title: str = "") -> Figure:
    """Plot time-domain analysis: rays, forces, torques.

    Reproduces Thesis Fig. 2-5 / 2-6 style plots.

    Parameters
    ----------
    results : dict
        Output from dynamics.run_time_domain().
    title : str
        Figure title.

    Returns
    -------
    fig : Figure
    """
    t = results["time"]
    q = results["q"]
    F = results["F_N"]
    N = q.shape[1]

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    # (a) Ray motion profile
    ax = axes[0]
    for i in range(N):
        ax.plot(t, np.rad2deg(q[:, i]), linewidth=1.5, label=f"$q_{{{i+1}}}$")
    if N >= 2:
        dq12 = np.rad2deg(q[:, 0] - q[:, 1])
        ax.plot(t, dq12, ":", color="gray", label=r"$q_1 - q_2$")
    ax.set_ylabel("Angle [°]")
    ax.set_title("(a) Motion profile of fin rays")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # (b) Hydrodynamic forces
    ax = axes[1]
    labels = [r"$F_{\mathrm{Heave}}$", r"$F_{\mathrm{Sway}}$", r"$F_{\mathrm{Surge}}$"]
    styles = ["-", "--", ":"]
    for k in range(3):
        ax.plot(t, F[:, k], styles[k], color="k", linewidth=1.5, label=labels[k])
    ax.set_ylabel("Force [N]")
    ax.set_title("(b) Hydrodynamic force components")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # (c) Hydrodynamic torques
    ax = axes[2]
    T_phi = results.get("T_phi")
    if T_phi is not None:
        for i in range(min(N, T_phi.shape[1])):
            ax.plot(t, T_phi[:, i], linewidth=1.5,
                    label=rf"$\tau_{{\Phi,{i+1}}}$")
        ax.set_ylabel("Torque [N·m]")
        ax.set_title("(c) Joint torques due to hydrodynamic loads")
        ax.legend(loc="upper right")
        ax.grid(True, alpha=0.3)
    ax.set_xlabel("Time [s]")

    if title:
        fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    return fig


def plot_parametric_1d(results: dict, xlabel: str = "", convert_deg: bool = False) -> Figure:
    """Plot 1-D parametric sweep results.

    Reproduces Thesis Fig. 5-3, 5-5 style.

    Parameters
    ----------
    results : dict
        Output from parametric.sweep_1d().
    xlabel : str
        Label for x-axis.
    convert_deg : bool
        Convert x-values from radians to degrees.
    """
    vals = results["values"]
    F = results["F_max"]
    pname = results["param_name"]

    if convert_deg:
        vals = np.rad2deg(vals)

    fig, ax = plt.subplots(figsize=(8, 5))
    labels = [r"$F_{\mathrm{Heave}}$", r"$F_{\mathrm{Sway}}$", r"$F_{\mathrm{Surge}}$"]
    styles = ["-", ":", "-."]
    for k in range(3):
        ax.plot(vals, F[:, k], styles[k], color="k", linewidth=1.5, label=labels[k])

    # Mark max surge
    idx_max = np.argmax(F[:, 2])
    ax.plot(vals[idx_max], F[idx_max, 2], "k*", markersize=10,
            label=f"max surge @ {vals[idx_max]:.1f}")

    ax.set_xlabel(xlabel or pname)
    ax.set_ylabel("Force [N]")
    ax.set_title(f"Parametric sweep: {pname}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def plot_ray_count_comparison(results: dict) -> Figure:
    """Plot forces vs N for different ray counts.

    Reproduces Thesis Fig. 5-15, 5-16 style.
    """
    phi = results["phi_values"]
    N_vals = results["N_values"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for title, data, ax in zip(
        ["Heave", "Sway", "Surge"],
        [results["F_heave"], results["F_sway"], results["F_surge"]],
        axes,
    ):
        for j, N in enumerate(N_vals):
            ax.plot(np.rad2deg(phi), data[:, j], linewidth=1.5, label=f"N={N}")
        ax.set_xlabel(r"$\phi$ [°]")
        ax.set_ylabel(f"$F_{{\\mathrm{{{title}}}}}$ [N]")
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle("Effect of number of rays on hydrodynamic forces", fontsize=14)
    fig.tight_layout()
    return fig


def plot_robot_motion(results: dict) -> Figure:
    """Plot robot position, velocity, acceleration.

    Reproduces Thesis Fig. 6-4 style.
    """
    t = results["time"]
    pos = results["pos"]
    vel = results["vel"]
    acc = results["acc"]

    fig, axes = plt.subplots(3, 3, figsize=(14, 10), sharex=True)
    labels = ["Surge (x)", "Sway (y)", "Heave (z)"]
    units = ["[m]", "[m/s]", "[m/s²]"]
    data = [pos, vel, acc]
    row_labels = ["Position", "Velocity", "Acceleration"]

    for row in range(3):
        for col in range(3):
            ax = axes[row, col]
            ax.plot(t, data[row][:, col], linewidth=1.5)
            ax.set_ylabel(f"{row_labels[row]} {units[row]}")
            ax.grid(True, alpha=0.3)
            if row == 0:
                ax.set_title(labels[col])
            if row == 2:
                ax.set_xlabel("Time [s]")

    fig.suptitle("Robot body dynamics", fontsize=14)
    fig.tight_layout()
    return fig


# ════════════════════════════════════════════════════════════════════
# Force stability / N-ray plots
# ════════════════════════════════════════════════════════════════════

def plot_N_stability(results: dict, title: str = "") -> Figure:
    """Plot force stability metrics vs number of rays N.

    Shows how increasing N improves force smoothness and reduces
    sign reversals (force oscillation ±).

    Parameters
    ----------
    results : dict
        Output from parametric.sweep_N_stability().
    """
    N_vals = results["N_values"]
    metrics = results["metrics"]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    comps = ["Heave", "Sway", "Surge"]
    colors = ["#2196F3", "#4CAF50", "#FF9800"]

    # Row 1: Ripple (%) and sign changes vs N
    for k, comp in enumerate(comps):
        ax = axes[0, k]
        ripples = [m[comp]["ripple_pct"] for m in metrics]
        sign_ch = [m[comp]["sign_changes"] for m in metrics]

        # Cap infinite ripple for display
        max_finite = max((r for r in ripples if r < float("inf")), default=100)
        ripples_disp = [r if r < float("inf") else max_finite * 1.5 for r in ripples]

        ax.plot(N_vals, ripples_disp, "o-", color=colors[k], lw=2, label="Ripple %")
        ax.set_ylabel("Ripple [%]", color=colors[k])
        ax.tick_params(axis="y", labelcolor=colors[k])

        ax2 = ax.twinx()
        ax2.plot(N_vals, sign_ch, "s--", color="gray", lw=1.5, label="Sign changes")
        ax2.set_ylabel("Sign changes", color="gray")
        ax2.tick_params(axis="y", labelcolor="gray")

        ax.set_xlabel("Number of rays $N$")
        ax.set_title(comp)
        ax.set_xticks(N_vals)
        ax.grid(True, alpha=0.3)

    # Row 2: Force smoothness (RMS dF/dt) and peak-to-peak
    for k, comp in enumerate(comps):
        ax = axes[1, k]
        smooth = [m[comp]["smoothness"] for m in metrics]
        p2p = [m[comp]["peak2peak"] for m in metrics]

        ax.plot(N_vals, smooth, "o-", color=colors[k], lw=2, label="RMS $dF/dt$")
        ax.set_ylabel("RMS $dF/dt$ [N/s]", color=colors[k])
        ax.tick_params(axis="y", labelcolor=colors[k])

        ax2 = ax.twinx()
        ax2.plot(N_vals, p2p, "s--", color="gray", lw=1.5, label="Peak-to-peak")
        ax2.set_ylabel("Peak-to-peak [N]", color="gray")
        ax2.tick_params(axis="y", labelcolor="gray")

        ax.set_xlabel("Number of rays $N$")
        ax.set_xticks(N_vals)
        ax.grid(True, alpha=0.3)

    fig.suptitle(title or "Force Stability vs Number of Rays $N$", fontsize=14)
    fig.tight_layout()
    return fig


def plot_N_force_waveforms(results: dict, t: np.ndarray,
                           N_highlight: list[int] | None = None) -> Figure:
    """Plot force time-series for selected N values overlaid.

    Visually shows how force waveforms become smoother with more rays.

    Parameters
    ----------
    results : dict
        Output from parametric.sweep_N_stability().
    t : array
        Time vector.
    N_highlight : list of int, optional
        Which N values to plot (default: all).
    """
    N_vals = results["N_values"]
    F_list = results["F_series"]
    comps = ["Heave", "Sway", "Surge"]

    if N_highlight is None:
        N_highlight = list(N_vals)

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(N_highlight)))

    for k, (comp, ax) in enumerate(zip(comps, axes)):
        for j, N in enumerate(N_vals):
            if int(N) not in N_highlight:
                continue
            cidx = N_highlight.index(int(N))
            F = F_list[j]
            ax.plot(t[:len(F)], F[:, k], color=cmap[cidx], lw=1.5,
                    alpha=0.8, label=f"$N={N}$")
        ax.set_ylabel(f"{comp} Force [N]")
        ax.legend(fontsize=8, loc="upper right", ncol=min(len(N_highlight), 4))
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time [s]")
    fig.suptitle("Force Waveforms: Effect of Ray Count $N$", fontsize=14)
    fig.tight_layout()
    return fig


def plot_geometry_stability(results: dict, param_label: str = "",
                            unit: str = "") -> Figure:
    """Plot stability metrics vs a swept physical dimension.

    Parameters
    ----------
    results : dict
        Output from parametric.sweep_geometry_stability().
    param_label : str
        Display label for x-axis.
    unit : str
        Unit string for x-axis.
    """
    vals = results["values"]
    metrics = results["metrics"]
    pname = results["param_name"]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    comps = ["Heave", "Sway", "Surge"]
    colors = ["#2196F3", "#4CAF50", "#FF9800"]

    for k, (comp, ax) in enumerate(zip(comps, axes)):
        smooth = [m[comp]["smoothness"] for m in metrics]
        p2p = [m[comp]["peak2peak"] for m in metrics]

        ax.plot(vals, smooth, "o-", color=colors[k], lw=2, label="RMS $dF/dt$")
        ax.set_ylabel("RMS $dF/dt$ [N/s]", color=colors[k])
        ax.tick_params(axis="y", labelcolor=colors[k])

        ax2 = ax.twinx()
        ax2.plot(vals, p2p, "s--", color="gray", lw=1.5, label="Peak-to-peak")
        ax2.set_ylabel("Peak-to-peak [N]", color="gray")
        ax2.tick_params(axis="y", labelcolor="gray")

        xlabel = param_label or pname
        if unit:
            xlabel += f" [{unit}]"
        ax.set_xlabel(xlabel)
        ax.set_title(comp)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"Force Stability vs {param_label or pname}", fontsize=14)
    fig.tight_layout()
    return fig
