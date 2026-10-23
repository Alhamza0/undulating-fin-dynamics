"""
Compare corrected Python analytical model against COMSOL CFD results.

Replicates the comparison from the MATLAB ``lastanlatcialvsCFD.m`` script,
but using the corrected Python implementation (fixed vn sign, elastic
torques) and with the MATLAB velocity bug fixed ([dqq1 dqq2] instead
of [dqq1 dqq1]).

CFD data files (in 3Rays_undFinDynamics/):
    F_phi0.xlsx  — Φ = 0°   (301×3, columns: [heave, sway, surge])
    F20.xlsx     — Φ = 20°  (301×3)
    F30.xlsx     — Φ = 30°  (343×3, trimmed to match time vector)
    F40.xlsx     — Φ = 40°  (301×3)
    F50.xlsx     — Φ = 50°  (301×3)

Column mapping (verified empirically against CFD):
    Python F[0] (sin terms) → Sway  (CFD col 1)
    Python F[1] (cos terms) → Heave (CFD col 0)
    Python F[2] (sin Δq)    → Surge (CFD col 2)

Parameters (matching COMSOL 2-ray model):
    A = 30°, f = 1 Hz, N_rays = 2, D = 0.09 m,
    hmin = 0.034 m, hmax = 0.174 m, Cn = 2.7, ρ = 1000 kg/m³
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Add the package to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from undulating_fin.params import FinParams
from undulating_fin.kinematics import compute_kinematics
from undulating_fin.hydrodynamics import compute_total_forces

# ─────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent.parent
CFD_DIR = BASE_DIR / "3Rays_undFinDynamics"
OUTPUT_DIR = Path(__file__).parent / "cfd_comparison_results"
OUTPUT_DIR.mkdir(exist_ok=True)

# Physical parameters (matching COMSOL 2-ray model / lastanlatcialvsCFD.m)
A_DEG = 30.0
F_HZ = 1.0
N_RAYS = 2
D = 0.09       # m
HMIN = 0.034   # m
HMAX = 0.174   # m
CN = 2.7
RHO = 1000.0   # kg/m³
X_DEG = 0.0    # angular offset

# Time parameters
N_CYCLES = 3
DT = 0.01 / F_HZ

# Phase shifts to compare (degrees)
PHI_CASES = {
    0:  "F_phi0.xlsx",
    20: "F20.xlsx",
    30: "F30.xlsx",
    40: "F40.xlsx",
    50: "F50.xlsx",
}

# Mapping from analytical Python indices to physical force components.
# Python compute_total_forces returns [F_x(sin), F_y(cos), F_z(sinΔq)].
# CFD (COMSOL) columns are [heave, sway, surge].
#   Physical Heave (vertical)   = Python F[1] (cos terms) → CFD col 0
#   Physical Sway  (lateral)    = Python F[0] (sin terms) → CFD col 1
#   Physical Surge (propulsive) = Python F[2] (sinΔq)     → CFD col 2
FORCE_MAP = [
    # (label, python_idx, cfd_idx)
    ("Heave", 1, 0),
    ("Sway",  0, 1),
    ("Surge", 2, 2),
]
FORCE_UNITS = "N"


def load_cfd_data(filename: str) -> np.ndarray:
    """Load CFD force data from Excel file.

    Files have no header; columns are [heave, sway, surge] in Newtons.
    Leading NaN rows (if any) are stripped.
    """
    filepath = CFD_DIR / filename
    df = pd.read_excel(filepath, header=None)
    data = df.values
    # Strip leading NaN rows (e.g. F30.xlsx has 42 NaN rows)
    valid_mask = ~np.any(np.isnan(data), axis=1)
    first_valid = np.argmax(valid_mask)
    return data[first_valid:]


def run_analytical(phi_deg: float, time_vec: np.ndarray) -> np.ndarray:
    """Run analytical model with prescribed sinusoidal kinematics.

    Returns raw array of shape (nT, 3) as from compute_total_forces:
    [F_x(sin), F_y(cos), F_z(sinΔq)].
    """
    params = FinParams(
        N=N_RAYS,
        A=np.deg2rad(A_DEG),
        f=F_HZ,
        Phi=np.deg2rad(phi_deg),
        X=np.deg2rad(X_DEG),
        D=D,
        hmax=HMAX,
        hmin=HMIN,
        Cn=CN,
        rho=RHO,
    )

    q, dq, _ = compute_kinematics(params, time_vec)

    nT = len(time_vec)
    F_analytical = np.zeros((nT, 3))
    for t_idx in range(nT):
        F_analytical[t_idx, :] = compute_total_forces(
            params, q[t_idx, :], dq[t_idx, :]
        )
        if (t_idx + 1) % 50 == 0 or t_idx == nT - 1:
            print(f"  Phi={phi_deg:3.0f}°: step {t_idx+1}/{nT}", end="\r")
    print()
    return F_analytical


def compute_metrics(
    analytical: np.ndarray, cfd: np.ndarray
) -> dict:
    """Compute comparison error metrics between analytical and CFD.

    Uses the last full cycle (steady-state) for metrics.
    Both inputs are already mapped: columns are [Heave, Sway, Surge].
    """
    steps_per_cycle = int(round(1.0 / (F_HZ * DT)))
    start_idx = max(0, len(analytical) - steps_per_cycle)
    an = analytical[start_idx:]
    cf = cfd[start_idx:]

    labels = [label for label, _, _ in FORCE_MAP]
    metrics = {}
    for k, label in enumerate(labels):
        diff = an[:, k] - cf[:, k]
        amp_cfd = np.max(np.abs(cf[:, k]))
        amp_cfd = max(amp_cfd, 1e-12)

        metrics[label] = {
            "max_abs_error": float(np.max(np.abs(diff))),
            "rmse": float(np.sqrt(np.mean(diff**2))),
            "nrmse_pct": float(np.sqrt(np.mean(diff**2)) / amp_cfd * 100),
            "correlation": float(np.corrcoef(an[:, k], cf[:, k])[0, 1])
                if np.std(cf[:, k]) > 1e-12 else float("nan"),
            "mean_an": float(np.mean(an[:, k])),
            "mean_cfd": float(np.mean(cf[:, k])),
        }
    return metrics


def remap_analytical(F_raw: np.ndarray) -> np.ndarray:
    """Remap raw analytical output to physical [Heave, Sway, Surge] order."""
    out = np.zeros_like(F_raw)
    for k, (_, py_idx, _) in enumerate(FORCE_MAP):
        out[:, k] = F_raw[:, py_idx]
    return out


def remap_cfd(cfd_raw: np.ndarray) -> np.ndarray:
    """Remap CFD columns to physical [Heave, Sway, Surge] order."""
    out = np.zeros_like(cfd_raw)
    for k, (_, _, cfd_idx) in enumerate(FORCE_MAP):
        out[:, k] = cfd_raw[:, cfd_idx]
    return out


def plot_comparison(
    time_vec: np.ndarray,
    analytical: np.ndarray,
    cfd: np.ndarray,
    phi_deg: float,
    metrics: dict,
):
    """Create a 3-subplot comparison figure (Heave, Sway, Surge).

    Both analytical and cfd are already in [Heave, Sway, Surge] order.
    """
    labels = [label for label, _, _ in FORCE_MAP]
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    fig.suptitle(
        f"Analytical vs CFD — Φ = {phi_deg}°\n"
        f"(N={N_RAYS} rays, A={A_DEG}°, f={F_HZ} Hz, D={D*100:.0f} cm, Cn={CN})",
        fontsize=13,
    )

    for k, (ax, label) in enumerate(zip(axes, labels)):
        ax.plot(time_vec, analytical[:, k], linewidth=2, label="Analytical (Python)")
        ax.plot(time_vec, cfd[:, k], "--", linewidth=2, label="CFD (COMSOL)")
        ax.set_ylabel(f"{label} [{FORCE_UNITS}]")
        ax.legend(loc="upper right", fontsize=9)
        ax.grid(True, alpha=0.3)

        m = metrics[label]
        info = f"RMSE={m['rmse']:.4f} N  NRMSE={m['nrmse_pct']:.1f}%  r={m['correlation']:.4f}"
        ax.text(
            0.02, 0.92, info,
            transform=ax.transAxes, fontsize=8,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8),
        )

    axes[-1].set_xlabel("Time [s]")
    plt.tight_layout()

    fname = OUTPUT_DIR / f"cfd_vs_analytical_phi{phi_deg:02.0f}.png"
    fig.savefig(fname, dpi=150)
    print(f"  Saved: {fname}")
    plt.close(fig)


def plot_all_phases_summary(all_results: dict):
    """Create summary figure comparing all phase shifts."""
    phi_list = sorted(all_results.keys())
    n_phi = len(phi_list)
    labels = [label for label, _, _ in FORCE_MAP]

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(
        f"Analytical vs CFD — All Phase Shifts\n"
        f"(N={N_RAYS} rays, A={A_DEG}°, f={F_HZ} Hz, D={D*100:.0f} cm, Cn={CN})",
        fontsize=13,
    )

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, n_phi))

    for k, (ax, label) in enumerate(zip(axes, labels)):
        for i, phi in enumerate(phi_list):
            r = all_results[phi]
            tv = r["time_vec"]
            ax.plot(tv, r["analytical"][:, k], color=colors[i],
                    linewidth=1.5, label=f"An. Φ={phi}°")
            ax.plot(tv, r["cfd"][:, k], "--", color=colors[i],
                    linewidth=1.5, alpha=0.7, label=f"CFD Φ={phi}°")
        ax.set_ylabel(f"{label} [{FORCE_UNITS}]")
        ax.grid(True, alpha=0.3)
        if k == 0:
            ax.legend(loc="upper right", fontsize=7, ncol=2)

    axes[-1].set_xlabel("Time [s]")
    plt.tight_layout()

    fname = OUTPUT_DIR / "cfd_vs_analytical_all_phases.png"
    fig.savefig(fname, dpi=150)
    print(f"  Saved: {fname}")
    plt.close(fig)


def plot_error_summary(all_results: dict):
    """Bar chart of NRMSE and correlation across phase shifts."""
    phi_list = sorted(all_results.keys())
    labels = [label for label, _, _ in FORCE_MAP]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Error Metrics: Analytical vs CFD (last cycle)", fontsize=13)

    x = np.arange(len(phi_list))
    width = 0.25

    for k, label in enumerate(labels):
        nrmse = [all_results[p]["metrics"][label]["nrmse_pct"] for p in phi_list]
        corr = [all_results[p]["metrics"][label]["correlation"] for p in phi_list]
        ax1.bar(x + k * width, nrmse, width, label=label)
        ax2.bar(x + k * width, corr, width, label=label)

    ax1.set_xticks(x + width)
    ax1.set_xticklabels([f"Φ={p}°" for p in phi_list])
    ax1.set_ylabel("NRMSE [%]")
    ax1.set_title("Normalized RMSE")
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis="y")

    ax2.set_xticks(x + width)
    ax2.set_xticklabels([f"Φ={p}°" for p in phi_list])
    ax2.set_ylabel("Correlation coefficient")
    ax2.set_title("Pearson Correlation")
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis="y")
    ax2.set_ylim(0.5, 1.02)

    plt.tight_layout()
    fname = OUTPUT_DIR / "cfd_error_metrics.png"
    fig.savefig(fname, dpi=150)
    print(f"  Saved: {fname}")
    plt.close(fig)


def main():
    print("=" * 70)
    print("  Analytical (Python) vs CFD (COMSOL) Comparison")
    print("=" * 70)
    print(f"  Parameters: A={A_DEG}°, f={F_HZ}Hz, N_rays={N_RAYS}, D={D}m")
    print(f"              hmin={HMIN}m, hmax={HMAX}m, Cn={CN}, rho={RHO}")
    print(f"  Phase shifts: {sorted(PHI_CASES.keys())}°")
    print(f"  Output dir: {OUTPUT_DIR}")
    print("=" * 70)

    # Time vector
    t_max = N_CYCLES / F_HZ
    time_vec = np.arange(0, t_max + DT / 2, DT)
    nT = len(time_vec)
    print(f"\n  Time steps: {nT}, dt={DT}s, tMax={t_max}s\n")

    all_results = {}
    labels = [label for label, _, _ in FORCE_MAP]

    for phi_deg, cfd_file in sorted(PHI_CASES.items()):
        print(f"\n{'─'*50}")
        print(f"  Processing Φ = {phi_deg}°  (CFD file: {cfd_file})")
        print(f"{'─'*50}")

        # Load CFD data and remap to [Heave, Sway, Surge]
        cfd_raw = load_cfd_data(cfd_file)
        print(f"  CFD data shape: {cfd_raw.shape}")

        # Trim to match analytical time vector length
        n_cfd = min(cfd_raw.shape[0], nT)
        cfd_mapped = remap_cfd(cfd_raw[:n_cfd, :])

        # Run analytical model and remap
        print(f"  Running analytical model...")
        F_raw = run_analytical(phi_deg, time_vec[:n_cfd])
        an_mapped = remap_analytical(F_raw)

        # Compute metrics on mapped data
        metrics = compute_metrics(an_mapped, cfd_mapped)

        # Print per-component metrics
        print(f"\n  {'Component':<10} {'RMSE [N]':>12} {'NRMSE [%]':>12} "
              f"{'Max Err [N]':>12} {'Correlation':>12}")
        print(f"  {'─'*58}")
        for label in labels:
            m = metrics[label]
            print(f"  {label:<10} {m['rmse']:>12.4f} {m['nrmse_pct']:>12.1f} "
                  f"{m['max_abs_error']:>12.4f} {m['correlation']:>12.4f}")

        all_results[phi_deg] = {
            "time_vec": time_vec[:n_cfd],
            "analytical": an_mapped,
            "cfd": cfd_mapped,
            "metrics": metrics,
        }

        plot_comparison(time_vec[:n_cfd], an_mapped, cfd_mapped, phi_deg, metrics)

    # Summary plots
    print(f"\n{'='*70}")
    print("  Generating summary plots...")
    plot_all_phases_summary(all_results)
    plot_error_summary(all_results)

    # Save metrics to CSV
    rows = []
    for phi_deg, r in sorted(all_results.items()):
        for label in labels:
            m = r["metrics"][label]
            rows.append({
                "Phi_deg": phi_deg,
                "Component": label,
                "RMSE_N": m["rmse"],
                "NRMSE_pct": m["nrmse_pct"],
                "Max_Abs_Error_N": m["max_abs_error"],
                "Correlation": m["correlation"],
                "Mean_Analytical_N": m["mean_an"],
                "Mean_CFD_N": m["mean_cfd"],
            })
    df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / "cfd_comparison_metrics.csv"
    df.to_csv(csv_path, index=False, float_format="%.6f")
    print(f"  Saved metrics: {csv_path}")

    print(f"\n{'='*70}")
    print("  Comparison complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
