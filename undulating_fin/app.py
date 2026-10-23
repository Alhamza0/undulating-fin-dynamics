"""
Undulating Fin Dynamics — Interactive Dashboard
================================================
Streamlit web application for kinematic and dynamic modeling
of biomimetic undulating fins for underwater robots.

Run with:  streamlit run app.py
"""

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from pathlib import Path
from scipy.signal import correlate

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from undulating_fin.params import FinParams, SimulationConfig, RobotParams
from undulating_fin.kinematics import compute_kinematics
from undulating_fin.dynamics import run_time_domain
from undulating_fin.hydrodynamics import compute_total_forces, compute_panel_forces
from undulating_fin.membrane import elastic_torques
from undulating_fin.parametric import sweep_1d, sweep_2d, sweep_ray_count, sweep_N_stability, sweep_geometry_stability, compute_stability_metrics
from undulating_fin.robot import MotionConfig, simulate_robot
from undulating_fin.plotting import (
    plot_time_domain,
    plot_parametric_1d,
    plot_ray_count_comparison,
    plot_robot_motion,
    plot_N_stability,
    plot_N_force_waveforms,
    plot_geometry_stability,
)

# OpenFOAM pipeline
from openfoam.case_generator import OFConfig, generate_case
from openfoam.result_parser import (
    load_forces, find_force_file, list_completed_cases,
    resample_to_analytical, ForceData,
)
from openfoam.validation import (
    validate_case, batch_validate, results_to_dataframe,
    run_analytical, FORCE_LABELS, compute_metrics,
)

# ──────────────────────────── page config ────────────────────────────
st.set_page_config(
    page_title="Undulating Fin Dynamics",
    page_icon="🐟",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ──────────────────────────── sidebar ────────────────────────────────
def build_sidebar():
    """Create sidebar with all fin parameters and return FinParams + SimulationConfig."""
    st.sidebar.title("⚙️ Parameters")

    # ── Wave / Kinematics ──
    with st.sidebar.expander("🌊 Wave & Kinematics", expanded=True):
        N = st.number_input("N — Number of rays", 2, 20, 3, step=1)
        A_deg = st.slider("A — Amplitude [°]", 1.0, 60.0, 30.0, 0.5)
        f = st.slider("f — Frequency [Hz]", 0.1, 6.0, 2.0, 0.1)
        Phi_deg = st.slider("Φ — Phase shift [°]", 0.0, 120.0, 40.0, 1.0)
        X_deg = st.slider("X — Bias angle [°]", -90.0, 90.0, 0.0, 1.0)

    # ── Geometry ──
    with st.sidebar.expander("📐 Geometry"):
        D_cm = st.slider("D — Ray spacing [cm]", 2.0, 20.0, 10.0, 0.5)
        hmax_cm = st.slider("h_max — Max membrane height [cm]", 1.0, 30.0, 17.5, 0.5)
        hmin_cm = st.slider("h_min — Min membrane height [cm]", 0.0, 20.0, 0.0, 0.5)

    # ── Fluid ──
    with st.sidebar.expander("💧 Fluid Properties"):
        rho = st.number_input("ρ — Fluid density [kg/m³]", 900.0, 1100.0, 1022.0, 1.0)
        Cn = st.slider("Cn — Normal drag coefficient", 1.0, 5.0, 2.65, 0.05)
        Cm = st.slider("Cm — Added mass coefficient (Morison)", 0.0, 3.0, 1.0, 0.05)

    # ── Hydrodynamic Model ──
    with st.sidebar.expander("🔬 Hydrodynamic Model", expanded=True):
        hydro_method = st.radio(
            "Force model",
            ["RFT (drag only)", "Morison (drag + inertia)", "Lighthill EBT (reactive)"],
            index=0,
            help=(
                "**RFT**: Resistive Force Theory — drag proportional to v_n². "
                "**Morison**: Adds added-mass inertia term (Cm). "
                "**Lighthill**: Large-amplitude elongated body theory — reactive forces."
            ),
        )
        method_map = {
            "RFT (drag only)": "rft",
            "Morison (drag + inertia)": "morison",
            "Lighthill EBT (reactive)": "lighthill",
        }
        hydro_key = method_map[hydro_method]

    # ── Membrane / Material ──
    with st.sidebar.expander("🔧 Membrane & Material"):
        EK = st.number_input("EK — Membrane stiffness [N/m]", 0.0, 1000.0, 360.0, 10.0)
        l_cm = st.slider("ℓ — Membrane char. length [cm]", 1.0, 20.0, 8.0, 0.5)
        v_poisson = st.slider("ν — Poisson's ratio", 0.0, 0.5, 0.5, 0.01)
        alpha = st.slider("α — Tension parameter D/ℓ", 0.1, 3.0, 1.0, 0.1)

    # ── Ray Mechanical ──
    with st.sidebar.expander("🦴 Ray Mechanical Properties"):
        mi = st.number_input("m_i — Ray mass [g]", 1.0, 100.0, 10.0, 1.0) / 1000.0
        li_cm = st.slider("l_i — Ray length [cm]", 3.0, 25.0, 12.0, 0.5)
        Izi = st.number_input("I_zi — Moment of inertia [×10⁻⁶ kg·m²]", 1.0, 500.0, 72.0, 1.0) * 1e-6
        bi = st.number_input("b_i — Viscous friction [N·m·s]", 0.0, 0.5, 0.02, 0.005, format="%.3f")
        kci = st.number_input("kc_i — Coulomb friction [N·m]", 0.0, 0.1, 0.001, 0.001, format="%.4f")

    # ── Simulation ──
    with st.sidebar.expander("⏱️ Simulation Settings"):
        t_max = st.slider("t_max — Duration [s]", 0.5, 20.0, 3.0, 0.5)
        dt = st.select_slider("dt — Time step [s]", options=[0.1, 0.05, 0.02, 0.01, 0.005], value=0.01)

    params = FinParams(
        N=int(N),
        A=np.deg2rad(A_deg),
        f=f,
        Phi=np.deg2rad(Phi_deg),
        X=np.deg2rad(X_deg),
        D=D_cm / 100.0,
        hmax=hmax_cm / 100.0,
        hmin=hmin_cm / 100.0,
        Cn=Cn,
        Cm=Cm,
        rho=rho,
        l=l_cm / 100.0,
        v=v_poisson,
        EK=EK,
        alpha=alpha,
        mi=mi,
        li=li_cm / 100.0,
        Izi=Izi,
        bi=bi,
        kci=kci,
    )
    sim = SimulationConfig(t_max=t_max, dt=dt)
    return params, sim, hydro_key


# ──────────────────────────── tabs ───────────────────────────────────

def tab_time_domain(params: FinParams, sim: SimulationConfig, hydro_method: str = "rft"):
    """Time-domain analysis — Thesis Sec. 2.4 / Fig. 2-5, 2-6."""
    st.header("Time-Domain Analysis")

    METHOD_LABELS = {"rft": "RFT (drag only)", "morison": "Morison (drag + inertia)", "lighthill": "Lighthill EBT (reactive)"}

    col1, col2 = st.columns([1.5, 1])
    with col2:
        st.markdown("#### Current Configuration")
        st.markdown(f"""
| Parameter | Value |
|-----------|-------|
| Rays (N) | **{params.N}** |
| Amplitude | **{np.rad2deg(params.A):.1f}°** |
| Frequency | **{params.f:.2f} Hz** |
| Phase shift | **{np.rad2deg(params.Phi):.1f}°** |
| Ray spacing | **{params.D*100:.1f} cm** |
| Wavelength | **{params.wavelength*100:.1f} cm** |
| Waves on fin | **{params.wave_count:.2f}** |
| Duration | **{sim.t_max:.1f} s** ({sim.n_steps} steps) |
| **Model** | **{METHOD_LABELS[hydro_method]}** |
""")

    with col1:
        compare_mode = st.checkbox("Compare all models", value=False, key="td_compare",
                                   help="Run RFT, Morison, and Lighthill side by side")
        if st.button("▶️ Run Simulation", key="run_td", type="primary", use_container_width=True):
            if compare_mode:
                all_results = {}
                for m in ("rft", "morison", "lighthill"):
                    with st.spinner(f"Running {METHOD_LABELS[m]}..."):
                        all_results[m] = run_time_domain(params, sim.time_vector, method=m)
                st.session_state["td_results"] = all_results["rft"]
                st.session_state["td_compare_results"] = all_results
            else:
                with st.spinner(f"Computing {METHOD_LABELS[hydro_method]} forces..."):
                    results = run_time_domain(params, sim.time_vector, method=hydro_method)
                st.session_state["td_results"] = results
                st.session_state.pop("td_compare_results", None)

    if "td_results" in st.session_state:
        results = st.session_state["td_results"]

        # ── Summary metrics ──
        F = results["F_N"]
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Max Heave Force", f"{np.max(np.abs(F[:,0])):.4f} N")
        c2.metric("Max Sway Force", f"{np.max(np.abs(F[:,1])):.4f} N")
        c3.metric("Max Surge Force", f"{np.max(np.abs(F[:,2])):.4f} N")
        c4.metric("Mean Surge Force", f"{np.mean(F[:,2]):.4f} N")

        # ── Plots ──
        fig = plot_time_domain(
            results,
            title=f"N={params.N}, A={np.rad2deg(params.A):.0f}°, "
                  f"Φ={np.rad2deg(params.Phi):.0f}°, f={params.f:.1f} Hz"
        )
        st.pyplot(fig)
        plt.close(fig)

        # ── Detailed torques ──
        with st.expander("📊 Detailed Torque Breakdown"):
            t = results["time"]
            fig2, axes = plt.subplots(1, 3, figsize=(15, 4))

            T_phi = results.get("T_phi")
            T_gamma = results.get("T_gamma")
            L = results.get("L")
            tau = results.get("tau")

            for i in range(min(params.N, 3)):
                if L is not None:
                    axes[0].plot(t, L[:, i], label=f"Ray {i+1}")
                if T_gamma is not None:
                    axes[1].plot(t, T_gamma[:, i], label=f"Ray {i+1}")
                if tau is not None:
                    axes[2].plot(t, tau[:, i], label=f"Ray {i+1}")

            axes[0].set_title("Inertia + Friction + Gravity (L)")
            axes[1].set_title("Elastic Membrane (Γ)")
            axes[2].set_title("Total Required Torque (τ)")
            for ax in axes:
                ax.set_xlabel("Time [s]")
                ax.set_ylabel("Torque [N·m]")
                ax.legend()
                ax.grid(True, alpha=0.3)
            fig2.tight_layout()
            st.pyplot(fig2)
            plt.close(fig2)

        # ── Raw data ──
        with st.expander("📋 Download Raw Data"):
            import pandas as pd
            df = pd.DataFrame({
                "Time [s]": results["time"],
                "F_heave [N]": F[:, 0],
                "F_sway [N]": F[:, 1],
                "F_surge [N]": F[:, 2],
            })
            # Add ray angles
            q = results["q"]
            for i in range(params.N):
                df[f"q{i+1} [deg]"] = np.rad2deg(q[:, i])
            st.dataframe(df, use_container_width=True, height=300)
            csv = df.to_csv(index=False)
            st.download_button("⬇️ Download CSV", csv, "time_domain_results.csv", "text/csv")

    # ── Model Comparison Plots ──
    if "td_compare_results" in st.session_state:
        st.divider()
        st.subheader("Model Comparison")
        cmp = st.session_state["td_compare_results"]
        t = cmp["rft"]["time"]

        METHOD_COLORS = {"rft": "#1f77b4", "morison": "#ff7f0e", "lighthill": "#2ca02c"}
        METHOD_LABELS_SHORT = {"rft": "RFT", "morison": "Morison", "lighthill": "Lighthill"}
        force_names = ["Heave (F_x)", "Sway (F_y)", "Surge (F_z)"]

        fig_cmp, axes = plt.subplots(1, 3, figsize=(15, 4))
        for k in range(3):
            for m in ("rft", "morison", "lighthill"):
                axes[k].plot(t, cmp[m]["F_N"][:, k],
                             color=METHOD_COLORS[m],
                             label=METHOD_LABELS_SHORT[m],
                             linewidth=1.5 if m == "rft" else 1.2)
            axes[k].set_title(force_names[k])
            axes[k].set_xlabel("Time [s]")
            axes[k].set_ylabel("Force [N]")
            axes[k].legend()
            axes[k].grid(True, alpha=0.3)
        fig_cmp.suptitle("Hydrodynamic Force Comparison: RFT vs Morison vs Lighthill", fontsize=12)
        fig_cmp.tight_layout()
        st.pyplot(fig_cmp)
        plt.close(fig_cmp)

        # ── Comparison metrics table ──
        rows = []
        for m in ("rft", "morison", "lighthill"):
            F_m = cmp[m]["F_N"]
            rows.append({
                "Model": METHOD_LABELS_SHORT[m],
                "Max |Heave| [N]": f"{np.max(np.abs(F_m[:, 0])):.5f}",
                "Max |Sway| [N]": f"{np.max(np.abs(F_m[:, 1])):.5f}",
                "Max |Surge| [N]": f"{np.max(np.abs(F_m[:, 2])):.5f}",
                "Mean Surge [N]": f"{np.mean(F_m[:, 2]):.5f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        # ── Difference from RFT ──
        with st.expander("Difference from RFT baseline"):
            F_rft = cmp["rft"]["F_N"]
            for m in ("morison", "lighthill"):
                F_m = cmp[m]["F_N"]
                diff = F_m - F_rft
                st.markdown(f"**{METHOD_LABELS_SHORT[m]} - RFT:**")
                col_a, col_b, col_c = st.columns(3)
                col_a.metric("Max Heave diff", f"{np.max(np.abs(diff[:, 0])):.5f} N")
                col_b.metric("Max Sway diff", f"{np.max(np.abs(diff[:, 1])):.5f} N")
                col_c.metric("Max Surge diff", f"{np.max(np.abs(diff[:, 2])):.5f} N")


def tab_parametric(params: FinParams, sim: SimulationConfig, hydro_method: str = "rft"):
    """Parametric sweeps — Thesis Ch. 5."""
    st.header("Parametric Sweep")
    st.caption("Sweep a single parameter while holding others fixed — reproduces Chapter 5 figures.")

    col_cfg, col_res = st.columns([1, 2])

    with col_cfg:
        sweep_type = st.radio(
            "Sweep type",
            ["Single Parameter (1D)", "Ray Count (N × Φ)", "N — Force Stability", "Geometry Stability"],
            horizontal=False,
        )

        if sweep_type == "Single Parameter (1D)":
            param_options = {
                "Φ — Phase shift": ("Phi", 0.0, 120.0, "°"),
                "A — Amplitude": ("A", 1.0, 60.0, "°"),
                "f — Frequency": ("f", 0.1, 6.0, "Hz"),
                "X — Bias angle": ("X", -90.0, 90.0, "°"),
                "D — Ray spacing": ("D", 2.0, 20.0, "cm"),
                "h_max — Membrane height": ("hmax", 1.0, 30.0, "cm"),
            }
            selected = st.selectbox("Parameter to sweep", list(param_options.keys()))
            pname, lo, hi, unit = param_options[selected]

            c1, c2 = st.columns(2)
            val_lo = c1.number_input(f"Start [{unit}]", value=lo, key="sw_lo")
            val_hi = c2.number_input(f"End [{unit}]", value=hi, key="sw_hi")
            n_pts = st.slider("Number of points", 5, 50, 15)

            if st.button("▶️ Run Sweep", key="run_1d", type="primary", use_container_width=True):
                raw_vals = np.linspace(val_lo, val_hi, n_pts)

                # Convert units for internal use
                if unit == "°":
                    internal_vals = np.deg2rad(raw_vals)
                elif unit == "cm":
                    internal_vals = raw_vals / 100.0
                else:
                    internal_vals = raw_vals

                with st.spinner(f"Sweeping {pname} ({n_pts} points)..."):
                    res = sweep_1d(params, sim, pname, internal_vals, verbose=False, method=hydro_method)
                st.session_state["sweep1d"] = res
                st.session_state["sweep1d_unit"] = unit
                st.session_state["sweep1d_label"] = selected

        elif sweep_type == "Ray Count (N × Φ)":
            c1, c2 = st.columns(2)
            n_min = c1.number_input("Min rays", 2, 12, 2)
            n_max = c2.number_input("Max rays", 2, 20, 6)
            n_phi = st.slider("Φ sweep points", 5, 40, 15, key="rc_pts")
            phi_lo = st.number_input("Φ start [°]", value=1.0, key="rc_lo")
            phi_hi = st.number_input("Φ end [°]", value=80.0, key="rc_hi")

            if st.button("▶️ Run Ray Count Sweep", key="run_rc", type="primary", use_container_width=True):
                N_vals = list(range(int(n_min), int(n_max) + 1))
                phi_vals = np.linspace(np.deg2rad(phi_lo), np.deg2rad(phi_hi), int(n_phi))
                with st.spinner(f"Sweeping N={N_vals} × {n_phi} Φ values..."):
                    res = sweep_ray_count(params, sim, N_vals, phi_vals, verbose=False, method=hydro_method)
                st.session_state["sweep_rc"] = res

        elif sweep_type == "N — Force Stability":
            st.markdown("""
            **How does the number of rays affect force smoothness?**

            More rays → wave closer to continuous sine → smoother forces
            → more stable robot motion.
            """)
            c1, c2 = st.columns(2)
            n_min = c1.number_input("Min N", 2, 6, 2, key="ns_min")
            n_max = c2.number_input("Max N", 4, 20, 12, key="ns_max")
            n_step = st.number_input("Step", 1, 5, 1, key="ns_step")

            if st.button("▶️ Run N-Stability Sweep", key="run_ns", type="primary", use_container_width=True):
                N_vals = list(range(int(n_min), int(n_max) + 1, int(n_step)))
                with st.spinner(f"Computing stability for N={N_vals}..."):
                    res = sweep_N_stability(params, sim, N_vals, verbose=False, method=hydro_method)
                st.session_state["sweep_ns"] = res

        else:  # Geometry Stability
            st.markdown("""
            **How do physical dimensions affect force stability?**

            Sweep ray spacing, membrane height, or ray length and
            analyze force smoothness metrics.
            """)
            geom_options = {
                "D — Ray spacing": ("D", 3.0, 20.0, "cm"),
                "h_max — Membrane height": ("hmax", 5.0, 30.0, "cm"),
                "l_i — Ray length": ("li", 5.0, 25.0, "cm"),
            }
            sel_geom = st.selectbox("Parameter", list(geom_options.keys()), key="gs_param")
            gpname, glo, ghi, gunit = geom_options[sel_geom]
            c1, c2 = st.columns(2)
            g_lo = c1.number_input(f"Start [{gunit}]", value=glo, key="gs_lo")
            g_hi = c2.number_input(f"End [{gunit}]", value=ghi, key="gs_hi")
            g_pts = st.slider("Points", 5, 30, 12, key="gs_pts")

            if st.button("▶️ Run Geometry Sweep", key="run_gs", type="primary", use_container_width=True):
                raw = np.linspace(g_lo, g_hi, g_pts)
                internal = raw / 100.0 if gunit == "cm" else raw
                with st.spinner(f"Sweeping {gpname} ({g_pts} points)..."):
                    res = sweep_geometry_stability(params, sim, gpname, internal,
                                                  verbose=False, method=hydro_method)
                st.session_state["sweep_gs"] = res
                st.session_state["sweep_gs_label"] = sel_geom.split("—")[0].strip()
                st.session_state["sweep_gs_unit"] = gunit

    with col_res:
        if "sweep1d" in st.session_state and sweep_type == "Single Parameter (1D)":
            res = st.session_state["sweep1d"]
            unit = st.session_state["sweep1d_unit"]
            lbl = st.session_state["sweep1d_label"]
            convert = unit == "°"
            fig = plot_parametric_1d(res, xlabel=lbl, convert_deg=convert)
            st.pyplot(fig)
            plt.close(fig)

            F = res["F_max"]
            idx_surge = np.argmax(F[:, 2])
            vals_display = np.rad2deg(res["values"]) if convert else res["values"]
            if unit == "cm":
                vals_display = res["values"] * 100

            st.success(f"**Optimal surge** at {res['param_name']} = {vals_display[idx_surge]:.1f} {unit}  →  "
                       f"F_surge = {F[idx_surge, 2]:.4f} N")

        if "sweep_rc" in st.session_state and sweep_type == "Ray Count (N × Φ)":
            res = st.session_state["sweep_rc"]
            fig = plot_ray_count_comparison(res)
            st.pyplot(fig)
            plt.close(fig)

        if "sweep_ns" in st.session_state and sweep_type == "N — Force Stability":
            res = st.session_state["sweep_ns"]

            # Stability metrics plot
            fig = plot_N_stability(res, title=f"Force Stability vs $N$ ({hydro_method.upper()})")
            st.pyplot(fig)
            plt.close(fig)

            # Waveform comparison
            fig2 = plot_N_force_waveforms(res, sim.time_vector)
            st.pyplot(fig2)
            plt.close(fig2)

            # Summary table
            st.subheader("Stability Summary")
            N_vals = res["N_values"]
            metrics = res["metrics"]
            rows = []
            for j, N in enumerate(N_vals):
                m = metrics[j]
                wn = res["wave_counts"][j]
                rows.append({
                    "N": int(N),
                    "Waves (wₙ)": f"{wn:.2f}",
                    "Surge sign changes": m["Surge"]["sign_changes"],
                    "Surge ripple %": f"{m['Surge']['ripple_pct']:.1f}" if m["Surge"]["ripple_pct"] < 1e6 else "∞",
                    "Surge RMS dF/dt": f"{m['Surge']['smoothness']:.4f}",
                    "Surge P2P [N]": f"{m['Surge']['peak2peak']:.6f}",
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        if "sweep_gs" in st.session_state and sweep_type == "Geometry Stability":
            res = st.session_state["sweep_gs"]
            lbl = st.session_state.get("sweep_gs_label", "")
            unt = st.session_state.get("sweep_gs_unit", "")
            fig = plot_geometry_stability(res, param_label=lbl, unit=unt)
            st.pyplot(fig)
            plt.close(fig)


def tab_robot(params: FinParams, sim: SimulationConfig):
    """Robot body dynamics — Thesis Ch. 6."""
    st.header("Robot Body Simulation")
    st.caption("Dual-fin underwater robot — reproduces Chapter 6 results.")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("🤖 Robot Body")
        R_cm = st.slider("R — Body radius [cm]", 3.0, 20.0, 10.0, 0.5)
        L_cm = st.slider("L — Body length [cm]", 10.0, 50.0, 22.0, 1.0)
        m_kg = st.number_input("m — Body mass [kg]", 0.5, 20.0, 2.5, 0.1)
        Cd_s = st.slider("Cd sphere", 0.5, 2.0, 1.1, 0.05)
        Cd_c = st.slider("Cd cylinder", 0.3, 1.5, 0.70, 0.05)

    with col2:
        st.subheader("🏊 Motion Mode")
        motion_type = st.radio("Select motion", ["Forward", "Backward", "Rotate CW", "Custom"],
                               horizontal=True)

        if motion_type == "Custom":
            st.markdown("**Right Fin**")
            r1, r2 = st.columns(2)
            rA = r1.slider("R — A [°]", 1.0, 60.0, 30.0, key="rA")
            rPhi = r2.slider("R — Φ [°]", 0.0, 120.0, 60.0, key="rPhi")

            st.markdown("**Left Fin**")
            l1, l2 = st.columns(2)
            lA = l1.slider("L — A [°]", 1.0, 60.0, 30.0, key="lA")
            lPhi = l2.slider("L — Φ [°]", -120.0, 120.0, -60.0, key="lPhi")
        else:
            rA = np.rad2deg(params.A)
            rPhi = np.rad2deg(params.Phi)

        sim_duration = st.slider("Simulation time [s]", 1.0, 30.0, 6.0, 0.5, key="robot_t")

    robot = RobotParams(
        R=R_cm / 100.0,
        L=L_cm / 100.0,
        m=m_kg,
        Cd_sphere=Cd_s,
        Cd_cylinder=Cd_c,
        rho=params.rho,
    )
    robot_sim = SimulationConfig(t_max=sim_duration, dt=sim.dt)

    if motion_type == "Forward":
        motion = MotionConfig.forward(A=params.A, f=params.f, Phi=params.Phi, N=params.N,
                                       D=params.D, hmax=params.hmax, hmin=params.hmin,
                                       Cn=params.Cn, rho=params.rho)
    elif motion_type == "Backward":
        motion = MotionConfig.backward(A=params.A, f=params.f, Phi=params.Phi, N=params.N,
                                        D=params.D, hmax=params.hmax, hmin=params.hmin,
                                        Cn=params.Cn, rho=params.rho)
    elif motion_type == "Rotate CW":
        motion = MotionConfig.rotate_cw(A=params.A, f=params.f, Phi=params.Phi, N=params.N,
                                         D=params.D, hmax=params.hmax, hmin=params.hmin,
                                         Cn=params.Cn, rho=params.rho)
    else:
        right_fin = FinParams(N=params.N, A=np.deg2rad(rA), f=params.f,
                              Phi=np.deg2rad(rPhi), X=params.X, D=params.D,
                              hmax=params.hmax, hmin=params.hmin,
                              Cn=params.Cn, rho=params.rho)
        left_fin = FinParams(N=params.N, A=np.deg2rad(lA), f=params.f,
                             Phi=np.deg2rad(lPhi), X=params.X, D=params.D,
                             hmax=params.hmax, hmin=params.hmin,
                             Cn=params.Cn, rho=params.rho)
        motion = MotionConfig(right_fin=right_fin, left_fin=left_fin,
                              description="Custom dual-fin configuration")

    if st.button("▶️ Simulate Robot", key="run_robot", type="primary", use_container_width=True):
        with st.spinner("Simulating robot body dynamics..."):
            results = simulate_robot(robot, motion, robot_sim)
        st.session_state["robot_results"] = results

    if "robot_results" in st.session_state:
        results = st.session_state["robot_results"]

        # Metrics
        pos = results["pos"]
        vel = results["vel"]
        c1, c2, c3 = st.columns(3)
        c1.metric("Final Surge (x)", f"{pos[-1, 0]*100:.2f} cm")
        c2.metric("Final Sway (y)", f"{pos[-1, 1]*100:.2f} cm")
        c3.metric("Final Heave (z)", f"{pos[-1, 2]*100:.2f} cm")

        c4, c5, c6 = st.columns(3)
        c4.metric("Peak Surge Vel.", f"{np.max(np.abs(vel[:,0])):.4f} m/s")
        c5.metric("Peak Sway Vel.", f"{np.max(np.abs(vel[:,1])):.4f} m/s")
        c6.metric("Peak Heave Vel.", f"{np.max(np.abs(vel[:,2])):.4f} m/s")

        # Plot
        fig = plot_robot_motion(results)
        st.pyplot(fig)
        plt.close(fig)

        # Force plots
        with st.expander("🔍 Fin Force Details"):
            t = results["time"]
            fig3, axes = plt.subplots(1, 3, figsize=(15, 4))
            labels = ["Heave", "Sway", "Surge"]
            for k in range(3):
                axes[k].plot(t, results["F_right"][:, k], label="Right fin", alpha=0.8)
                axes[k].plot(t, results["F_left"][:, k], label="Left fin", alpha=0.8)
                axes[k].plot(t, results["F_total"][:, k], "k--", label="Total", linewidth=2)
                axes[k].set_title(labels[k])
                axes[k].set_xlabel("Time [s]")
                axes[k].set_ylabel("Force [N]")
                axes[k].legend()
                axes[k].grid(True, alpha=0.3)
            fig3.tight_layout()
            st.pyplot(fig3)
            plt.close(fig3)


def tab_kinematics(params: FinParams, sim: SimulationConfig):
    """Interactive kinematics visualization."""
    st.header("Kinematics Visualization")
    st.caption("Visualize the traveling wave — ray angles, velocities, and accelerations (Eq. 1-3).")

    q, dq, ddq = compute_kinematics(params, sim.time_vector)
    t = sim.time_vector

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)

    # Angles
    for i in range(params.N):
        axes[0].plot(t, np.rad2deg(q[:, i]), linewidth=1.5, label=f"$q_{{{i+1}}}$")
    axes[0].set_ylabel("Angle [°]")
    axes[0].set_title("(a) Joint Angles — q(t)")
    axes[0].legend(ncol=min(params.N, 6), loc="upper right", fontsize=8)
    axes[0].grid(True, alpha=0.3)

    # Velocities
    for i in range(params.N):
        axes[1].plot(t, np.rad2deg(dq[:, i]), linewidth=1.5, label=f"q̇{i+1}")
    axes[1].set_ylabel("Velocity [°/s]")
    axes[1].set_title("(b) Joint Velocities — q̇(t)")
    axes[1].legend(ncol=min(params.N, 6), loc="upper right", fontsize=8)
    axes[1].grid(True, alpha=0.3)

    # Accelerations
    for i in range(params.N):
        axes[2].plot(t, np.rad2deg(ddq[:, i]), linewidth=1.5, label=f"q̈{i+1}")
    axes[2].set_ylabel("Acceleration [°/s²]")
    axes[2].set_title("(c) Joint Accelerations — q̈(t)")
    axes[2].set_xlabel("Time [s]")
    axes[2].legend(ncol=min(params.N, 6), loc="upper right", fontsize=8)
    axes[2].grid(True, alpha=0.3)

    fig.suptitle(
        f"Traveling Wave: N={params.N}, A={np.rad2deg(params.A):.0f}°, "
        f"f={params.f:.1f} Hz, Φ={np.rad2deg(params.Phi):.0f}°",
        fontsize=14
    )
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)

    # Wave info
    st.markdown(f"""
**Wave Properties:**
- Wavelength: λ = 2πD / Φ = **{params.wavelength*100:.2f} cm**
- Wave count on fin: **{params.wave_count:.2f}**
- Period: T = 1/f = **{1/params.f:.3f} s**
- Max angle: A + X = **{np.rad2deg(params.A + params.X):.1f}°**
""")


# ──────────────────── CFD Comparison Constants ───────────────────────
_CFD_DIR = Path(__file__).resolve().parent.parent / "3Rays_undFinDynamics"

_CFD_A_DEG = 30.0
_CFD_F_HZ = 1.0
_CFD_N_RAYS = 2
_CFD_D = 0.09
_CFD_HMIN = 0.034
_CFD_HMAX = 0.174
_CFD_CN = 2.7
_CFD_RHO = 1000.0
_CFD_DT = 0.01
_CFD_N_CYCLES = 3

# Bump this whenever the force computation code changes, to invalidate
# Streamlit's @st.cache_data results that depend on compute_total_forces.
_MODEL_VERSION = 2  # v2: fixed Morison inertia sign (body-in-still-fluid)

_CFD_PHI_CASES = {
    0: "F_phi0.xlsx", 20: "F20.xlsx", 30: "F30.xlsx",
    40: "F40.xlsx", 50: "F50.xlsx",
}

_FORCE_MAP = [("Heave", 1, 0), ("Sway", 0, 1), ("Surge", 2, 2)]

# Column indices for COMSOL kinematics files (qdq_step.xlsx / dq.xlsx)
_PHI_COLS = {
    20: (0, 1, 0, 1), 30: (4, 5, 2, 3),
    40: (8, 9, 4, 5), 50: (12, 13, 6, 7),
}


@st.cache_data(show_spinner=False)
def _load_cfd(filename: str) -> np.ndarray:
    """Load CFD force data; strip leading NaN rows."""
    data = pd.read_excel(_CFD_DIR / filename, header=None).values
    valid = ~np.any(np.isnan(data), axis=1)
    return data[np.argmax(valid):]


@st.cache_data(show_spinner=False)
def _run_cfd_analytical(phi_deg: float, *, method: str = "rft", Cm: float = 1.0,
                        _ver: int = _MODEL_VERSION):
    """Run analytical model and return remapped (Heave, Sway, Surge) array."""
    time_vec = np.arange(0, _CFD_N_CYCLES / _CFD_F_HZ + _CFD_DT / 2, _CFD_DT)
    params = FinParams(
        N=_CFD_N_RAYS, A=np.deg2rad(_CFD_A_DEG), f=_CFD_F_HZ,
        Phi=np.deg2rad(phi_deg), X=0.0, D=_CFD_D,
        hmax=_CFD_HMAX, hmin=_CFD_HMIN, Cn=_CFD_CN, Cm=Cm, rho=_CFD_RHO,
    )
    q, dq, ddq = compute_kinematics(params, time_vec)
    F_raw = np.zeros((len(time_vec), 3))
    for i in range(len(time_vec)):
        F_raw[i, :] = compute_total_forces(
            params, q[i, :], dq[i, :], ddq=ddq[i, :], method=method,
        )
    # Remap to physical order
    F = np.column_stack([F_raw[:, py] for _, py, _ in _FORCE_MAP])
    return time_vec, F


@st.cache_data(show_spinner=False)
def _tune_cm(phi_deg: float, cm_values: tuple[float, ...],
             _ver: int = _MODEL_VERSION):
    """Sweep Cm values for Morison model and compute NRMSE against CFD."""
    cfd_raw = _load_cfd(_CFD_PHI_CASES[phi_deg])
    F_cfd = np.column_stack([cfd_raw[:, ci] for _, _, ci in _FORCE_MAP])
    time_vec = np.arange(0, _CFD_N_CYCLES / _CFD_F_HZ + _CFD_DT / 2, _CFD_DT)
    n = min(len(time_vec), len(F_cfd))
    spc = int(round(1.0 / (_CFD_F_HZ * _CFD_DT)))
    s = max(0, n - spc)  # last cycle only

    results = {"cm": list(cm_values)}
    for lbl, _, _ in _FORCE_MAP:
        results[lbl] = {"nrmse": [], "corr": []}
    results["mean_nrmse"] = []

    for cm in cm_values:
        _, F_an = _run_cfd_analytical(phi_deg, method="morison", Cm=cm)
        F_an = F_an[:n]
        nrmses = []
        for k, (lbl, _, _) in enumerate(_FORCE_MAP):
            ref = F_cfd[s:n, k]
            sig = F_an[s:n, k]
            amp = max(np.max(np.abs(ref)), 1e-12)
            nrmse = float(np.sqrt(np.mean((sig - ref) ** 2)) / amp * 100)
            corr = float(np.corrcoef(sig, ref)[0, 1]) if np.std(ref) > 1e-12 else float("nan")
            results[lbl]["nrmse"].append(nrmse)
            results[lbl]["corr"].append(corr)
            nrmses.append(nrmse)
        results["mean_nrmse"].append(float(np.mean(nrmses)))

    return results


@st.cache_data(show_spinner=False)
def _run_phase_lag_analysis():
    """Compute forces with [dq1,dq2] vs [dq1,dq1] for all phase shifts."""
    qdq = pd.read_excel(_CFD_DIR / "qdq_step.xlsx", header=None).values
    dq_all = pd.read_excel(_CFD_DIR / "dq.xlsx", header=None).values
    time_vec = np.arange(0, _CFD_N_CYCLES / _CFD_F_HZ + _CFD_DT / 2, _CFD_DT)
    nT = len(time_vec)

    results = {}
    for phi_deg in [20, 30, 40, 50]:
        cfd = _load_cfd(_CFD_PHI_CASES[phi_deg])
        qc1, qc2, dc1, dc2 = _PHI_COLS[phi_deg]
        n = min(nT, len(cfd))

        F_correct = np.zeros((n, 3))
        F_match = np.zeros((n, 3))
        for i in range(n):
            q1, q2 = qdq[i, qc1], qdq[i, qc2]
            d1, d2 = dq_all[i, dc1], dq_all[i, dc2]
            F_correct[i, :] = compute_panel_forces(
                _CFD_D, _CFD_HMIN, _CFD_HMAX, q1, q2, d1, d2, _CFD_CN, _CFD_RHO)
            F_match[i, :] = compute_panel_forces(
                _CFD_D, _CFD_HMIN, _CFD_HMAX, q1, q2, d1, d1, _CFD_CN, _CFD_RHO)

        skip = 100  # skip first cycle transient
        phi_metrics = {}
        for lbl, f_idx, c_idx in _FORCE_MAP:
            cf = cfd[skip:n, c_idx]
            an_c = F_correct[skip:n, f_idx]
            an_m = F_match[skip:n, f_idx]

            def _metrics(sig, ref):
                r = float(np.corrcoef(sig, ref)[0, 1]) if np.std(ref) > 1e-12 else float('nan')
                s_n = (sig - np.mean(sig)) / (np.std(sig) + 1e-12)
                r_n = (ref - np.mean(ref)) / (np.std(ref) + 1e-12)
                xc = correlate(r_n, s_n, mode='full')
                lags = np.arange(-len(s_n) + 1, len(s_n))
                lag = float(lags[np.argmax(xc)] * _CFD_DT)
                nrmse = float(np.sqrt(np.mean((sig - ref) ** 2))
                              / max(np.max(np.abs(ref)), 1e-12) * 100)
                return r, lag, nrmse

            r1, l1, n1 = _metrics(an_c, cf)
            r2, l2, n2 = _metrics(an_m, cf)
            phi_metrics[lbl] = {
                'correct': {'corr': r1, 'lag': l1, 'nrmse': n1},
                'match': {'corr': r2, 'lag': l2, 'nrmse': n2},
            }

        results[phi_deg] = {
            'F_correct': F_correct, 'F_match': F_match,
            'cfd': cfd[:n], 'n': n, 'metrics': phi_metrics,
        }
    return results, time_vec


def tab_cfd_comparison():
    """CFD Comparison — analytical model vs COMSOL FSI results."""
    st.header("CFD Comparison & Validation")
    st.caption("Compare the analytical model with COMSOL FSI simulation results across multiple phase shifts.")

    # ─── Check data availability ───
    if not _CFD_DIR.exists():
        st.error(f"CFD data directory not found: `{_CFD_DIR}`")
        return

    # ─── Section 1: Methodology ───
    with st.expander("📖 Methodology & Parameters", expanded=False):
        st.markdown("""
#### Model Validation Approach

The analytical model (Resistive Force Theory — RFT) is validated against high-fidelity
**COMSOL Multiphysics** Fluid-Structure Interaction (FSI) simulations. The COMSOL model
solves the full Navier–Stokes equations coupled with a moving mesh for the undulating fin
geometry.

**Comparison Setup:**
- **Geometry:** 2-ray fin (N=2) with inter-ray spacing D = 9 cm
- **Kinematics:** Sinusoidal oscillation, A = 30°, f = 1 Hz
- **Membrane:** hmin = 3.4 cm, hmax = 17.4 cm
- **Fluid:** Water (ρ = 1000 kg/m³), drag coefficient Cn = 2.7
- **Phase shifts tested:** Φ = 0°, 20°, 30°, 40°, 50°
- **Simulation duration:** 3 seconds (3 cycles) at dt = 0.01 s
- **Metrics computed on final (steady-state) cycle**

**Force Components:**
| Component | Physical Meaning | Python Index | CFD Column |
|-----------|-----------------|:------------:|:----------:|
| Heave | Vertical (lift) force | F[1] (cos terms) | Column 0 |
| Sway | Lateral force | F[0] (sin terms) | Column 1 |
| Surge | Propulsive (thrust) force | F[2] (sin Δq terms) | Column 2 |
""")

    # ─── Section 2: Interactive Comparison ───
    st.subheader("Force Comparison")

    col_ctrl, col_plot = st.columns([1, 3])

    with col_ctrl:
        phi_select = st.selectbox(
            "Phase shift Φ [°]", list(_CFD_PHI_CASES.keys()),
            index=3, key="cfd_phi_select",
        )
        show_all = st.checkbox("Show all phase shifts", value=False, key="cfd_show_all")
        run_cfd = st.button("▶️ Run Comparison", key="run_cfd_btn",
                            type="primary", use_container_width=True)

    if run_cfd:
        if show_all:
            all_results = {}
            with st.spinner("Computing analytical forces for all phase shifts..."):
                for phi in _CFD_PHI_CASES:
                    tv, F_an = _run_cfd_analytical(phi)
                    cfd_raw = _load_cfd(_CFD_PHI_CASES[phi])
                    F_cfd = np.column_stack([cfd_raw[:, ci] for _, _, ci in _FORCE_MAP])
                    n = min(len(tv), len(F_cfd))
                    spc = int(round(1.0 / (_CFD_F_HZ * _CFD_DT)))
                    s = max(0, n - spc)
                    metrics = {}
                    for k, (lbl, _, _) in enumerate(_FORCE_MAP):
                        diff = F_an[s:n, k] - F_cfd[s:n, k]
                        amp = max(np.max(np.abs(F_cfd[s:n, k])), 1e-12)
                        rmse = float(np.sqrt(np.mean(diff ** 2)))
                        r_val = float(np.corrcoef(F_an[s:n, k], F_cfd[s:n, k])[0, 1]) \
                            if np.std(F_cfd[s:n, k]) > 1e-12 else float('nan')
                        metrics[lbl] = {'rmse': rmse, 'nrmse': rmse / amp * 100, 'corr': r_val}
                    all_results[phi] = {
                        'time': tv[:n], 'analytical': F_an[:n], 'cfd': F_cfd[:n], 'metrics': metrics
                    }
            st.session_state["cfd_all"] = all_results
        else:
            with st.spinner(f"Computing analytical forces for Φ = {phi_select}°..."):
                tv, F_an = _run_cfd_analytical(phi_select)
                cfd_raw = _load_cfd(_CFD_PHI_CASES[phi_select])
                F_cfd = np.column_stack([cfd_raw[:, ci] for _, _, ci in _FORCE_MAP])
                n = min(len(tv), len(F_cfd))
                spc = int(round(1.0 / (_CFD_F_HZ * _CFD_DT)))
                s = max(0, n - spc)
                metrics = {}
                for k, (lbl, _, _) in enumerate(_FORCE_MAP):
                    diff = F_an[s:n, k] - F_cfd[s:n, k]
                    amp = max(np.max(np.abs(F_cfd[s:n, k])), 1e-12)
                    rmse = float(np.sqrt(np.mean(diff ** 2)))
                    r_val = float(np.corrcoef(F_an[s:n, k], F_cfd[s:n, k])[0, 1]) \
                        if np.std(F_cfd[s:n, k]) > 1e-12 else float('nan')
                    metrics[lbl] = {'rmse': rmse, 'nrmse': rmse / amp * 100, 'corr': r_val}
            st.session_state["cfd_single"] = {
                'phi': phi_select, 'time': tv[:n],
                'analytical': F_an[:n], 'cfd': F_cfd[:n], 'metrics': metrics
            }

    with col_plot:
        # Single phase shift plot
        if "cfd_single" in st.session_state and not show_all:
            r = st.session_state["cfd_single"]
            fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
            fig.suptitle(
                f"Analytical vs CFD — Φ = {r['phi']}°\n"
                f"(N=2, A=30°, f=1 Hz, D=9 cm, Cn=2.7)", fontsize=12)
            labels = [lbl for lbl, _, _ in _FORCE_MAP]
            for k, (ax, lbl) in enumerate(zip(axes, labels)):
                ax.plot(r['time'], r['analytical'][:, k], lw=2, label='Analytical (RFT)')
                ax.plot(r['time'], r['cfd'][:, k], '--', lw=2, label='CFD (COMSOL)')
                ax.set_ylabel(f"{lbl} [N]")
                ax.legend(loc='upper right', fontsize=9)
                ax.grid(True, alpha=0.3)
                m = r['metrics'][lbl]
                ax.text(0.02, 0.92,
                        f"RMSE={m['rmse']:.4f} N  NRMSE={m['nrmse']:.1f}%  r={m['corr']:.3f}",
                        transform=ax.transAxes, fontsize=8, va='top',
                        bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.8))
            axes[-1].set_xlabel("Time [s]")
            fig.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

        # All phases plot
        if "cfd_all" in st.session_state and show_all:
            all_r = st.session_state["cfd_all"]
            phi_list = sorted(all_r.keys())
            colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(phi_list)))
            labels = [lbl for lbl, _, _ in _FORCE_MAP]

            fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
            fig.suptitle("Analytical vs CFD — All Phase Shifts", fontsize=12)
            for k, (ax, lbl) in enumerate(zip(axes, labels)):
                for i, phi in enumerate(phi_list):
                    d = all_r[phi]
                    ax.plot(d['time'], d['analytical'][:, k], color=colors[i],
                            lw=1.5, label=f'An. Φ={phi}°')
                    ax.plot(d['time'], d['cfd'][:, k], '--', color=colors[i],
                            lw=1.5, alpha=0.7, label=f'CFD Φ={phi}°')
                ax.set_ylabel(f"{lbl} [N]")
                ax.grid(True, alpha=0.3)
                if k == 0:
                    ax.legend(loc='upper right', fontsize=7, ncol=2)
            axes[-1].set_xlabel("Time [s]")
            fig.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

    # ─── Metrics Table ───
    if "cfd_all" in st.session_state and show_all:
        st.subheader("Error Metrics (Last Cycle)")
        all_r = st.session_state["cfd_all"]
        rows = []
        for phi in sorted(all_r.keys()):
            for lbl in ["Heave", "Sway", "Surge"]:
                m = all_r[phi]['metrics'][lbl]
                rows.append({
                    "Φ [°]": phi, "Component": lbl,
                    "RMSE [N]": f"{m['rmse']:.4f}",
                    "NRMSE [%]": f"{m['nrmse']:.1f}",
                    "Correlation": f"{m['corr']:.3f}",
                })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

        # Error bar chart
        phi_list = sorted(all_r.keys())
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle("Error Metrics Summary (Last Cycle)", fontsize=12)
        x = np.arange(len(phi_list))
        w = 0.25
        for k, lbl in enumerate(["Heave", "Sway", "Surge"]):
            nrmse = [all_r[p]['metrics'][lbl]['nrmse'] for p in phi_list]
            corr = [all_r[p]['metrics'][lbl]['corr'] for p in phi_list]
            ax1.bar(x + k * w, nrmse, w, label=lbl)
            ax2.bar(x + k * w, corr, w, label=lbl)
        ax1.set_xticks(x + w)
        ax1.set_xticklabels([f"Φ={p}°" for p in phi_list])
        ax1.set_ylabel("NRMSE [%]")
        ax1.set_title("Normalized RMSE")
        ax1.legend()
        ax1.grid(True, alpha=0.3, axis='y')
        ax2.set_xticks(x + w)
        ax2.set_xticklabels([f"Φ={p}°" for p in phi_list])
        ax2.set_ylabel("Correlation")
        ax2.set_title("Pearson Correlation")
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.set_ylim(0.5, 1.02)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)

    elif "cfd_single" in st.session_state and not show_all:
        st.subheader("Error Metrics (Last Cycle)")
        r = st.session_state["cfd_single"]
        rows = []
        for lbl in ["Heave", "Sway", "Surge"]:
            m = r['metrics'][lbl]
            rows.append({
                "Component": lbl,
                "RMSE [N]": f"{m['rmse']:.4f}",
                "NRMSE [%]": f"{m['nrmse']:.1f}",
                "Correlation": f"{m['corr']:.3f}",
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # ─── Section 3: Cm Auto-Tuning ───
    st.divider()
    st.subheader("Morison Cm Tuning — Match CFD")
    st.caption(
        "Sweep the added-mass coefficient Cm to find the value that best "
        "matches CFD forces. The optimal Cm minimizes the mean NRMSE across "
        "all three force components on the last (steady-state) cycle."
    )

    col_t1, col_t2 = st.columns([1, 1])
    with col_t1:
        tune_phi = st.selectbox(
            "Phase shift for tuning",
            list(_CFD_PHI_CASES.keys()),
            index=3,
            key="tune_phi",
        )
    with col_t2:
        cm_range = st.slider(
            "Cm search range", 0.0, 5.0, (0.0, 3.0), 0.1, key="cm_range",
        )
    cm_steps = st.slider("Number of Cm samples", 10, 60, 25, key="cm_steps")

    if st.button("🔍 Find Optimal Cm", key="tune_cm_btn", type="primary"):
        cm_vals = tuple(np.linspace(cm_range[0], cm_range[1], cm_steps).tolist())
        with st.spinner(f"Sweeping {cm_steps} Cm values for Φ = {tune_phi}°..."):
            tune_res = _tune_cm(tune_phi, cm_vals)
        st.session_state["tune_res"] = tune_res
        st.session_state["tune_phi_val"] = tune_phi

    if "tune_res" in st.session_state:
        tr = st.session_state["tune_res"]
        cm_arr = np.array(tr["cm"])
        mean_nrmse = np.array(tr["mean_nrmse"])
        best_idx = int(np.argmin(mean_nrmse))
        best_cm = cm_arr[best_idx]
        best_err = mean_nrmse[best_idx]
        tune_phi_val = st.session_state.get("tune_phi_val", "?")

        st.success(f"**Optimal Cm = {best_cm:.2f}** — Mean NRMSE = {best_err:.1f}%  (Φ = {tune_phi_val}°)")

        # ── NRMSE vs Cm plot ──
        fig_t, ax_t = plt.subplots(figsize=(10, 4))
        ax_t.plot(cm_arr, mean_nrmse, "k-o", ms=4, lw=2, label="Mean NRMSE")
        colors_comp = {"Heave": "#1f77b4", "Sway": "#ff7f0e", "Surge": "#2ca02c"}
        for lbl, _, _ in _FORCE_MAP:
            ax_t.plot(cm_arr, tr[lbl]["nrmse"], "--", lw=1.2,
                      color=colors_comp[lbl], alpha=0.7, label=lbl)
        ax_t.axvline(best_cm, color="red", ls=":", lw=1.5, label=f"Best Cm = {best_cm:.2f}")
        ax_t.set_xlabel("Cm (Added Mass Coefficient)")
        ax_t.set_ylabel("NRMSE [%]")
        ax_t.set_title(f"Morison Model Error vs Cm  (Φ = {tune_phi_val}°, last cycle)")
        ax_t.legend(fontsize=9)
        ax_t.grid(True, alpha=0.3)
        fig_t.tight_layout()
        st.pyplot(fig_t)
        plt.close(fig_t)

        # ── Correlation vs Cm plot ──
        fig_c, ax_c = plt.subplots(figsize=(10, 4))
        for lbl, _, _ in _FORCE_MAP:
            ax_c.plot(cm_arr, tr[lbl]["corr"], "-s", ms=3, lw=1.5,
                      color=colors_comp[lbl], label=lbl)
        ax_c.axvline(best_cm, color="red", ls=":", lw=1.5)
        ax_c.set_xlabel("Cm")
        ax_c.set_ylabel("Pearson Correlation")
        ax_c.set_title("Correlation with CFD vs Cm")
        ax_c.legend(fontsize=9)
        ax_c.grid(True, alpha=0.3)
        ax_c.set_ylim(0, 1.05)
        fig_c.tight_layout()
        st.pyplot(fig_c)
        plt.close(fig_c)

        # ── Best vs RFT vs CFD overlay ──
        with st.expander("📊 Overlay: RFT vs Morison (best Cm) vs CFD"):
            _, F_rft = _run_cfd_analytical(tune_phi_val, method="rft")
            _, F_mor = _run_cfd_analytical(tune_phi_val, method="morison", Cm=best_cm)
            cfd_raw = _load_cfd(_CFD_PHI_CASES[tune_phi_val])
            F_cfd = np.column_stack([cfd_raw[:, ci] for _, _, ci in _FORCE_MAP])
            tv = np.arange(0, _CFD_N_CYCLES / _CFD_F_HZ + _CFD_DT / 2, _CFD_DT)
            n = min(len(tv), len(F_cfd), len(F_rft), len(F_mor))

            fig_ov, axes_ov = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
            fig_ov.suptitle(
                f"RFT vs Morison (Cm={best_cm:.2f}) vs CFD — Φ = {tune_phi_val}°",
                fontsize=12,
            )
            labels_ov = [lbl for lbl, _, _ in _FORCE_MAP]
            for k, (ax, lbl) in enumerate(zip(axes_ov, labels_ov)):
                ax.plot(tv[:n], F_cfd[:n, k], "k-", lw=2, label="CFD (COMSOL)")
                ax.plot(tv[:n], F_rft[:n, k], "--", lw=1.5, color="#1f77b4", label="RFT")
                ax.plot(tv[:n], F_mor[:n, k], "-.", lw=1.5, color="#ff7f0e",
                        label=f"Morison (Cm={best_cm:.2f})")
                ax.set_ylabel(f"{lbl} [N]")
                ax.legend(loc="upper right", fontsize=9)
                ax.grid(True, alpha=0.3)
            axes_ov[-1].set_xlabel("Time [s]")
            fig_ov.tight_layout()
            st.pyplot(fig_ov)
            plt.close(fig_ov)

        st.info(
            f"💡 **Recommendation:** Set Cm = **{best_cm:.2f}** in the sidebar "
            f"(Fluid Properties → Cm slider) to use this tuned value in other tabs."
        )

    # ─── Section 4: Phase Lag / Velocity Analysis ───
    st.divider()
    st.subheader("Phase Lag Analysis — COMSOL Velocity Limitation")

    st.markdown("""
The primary source of discrepancy between analytical and CFD results was traced to the
**COMSOL FSI model's moving mesh boundary condition**. Specifically, the moving mesh
derived membrane velocity exclusively from **joint 1** ($\\dot{q}_1$), ignoring the
independent velocity contribution of **joint 2** ($\\dot{q}_2$).

This means the COMSOL simulation effectively used $[\\dot{q}_1, \\dot{q}_1]$ instead of
the correct $[\\dot{q}_1, \\dot{q}_2]$. To verify this, forces are computed both ways
and compared with the CFD output.
""")

    if not (_CFD_DIR / "qdq_step.xlsx").exists() or not (_CFD_DIR / "dq.xlsx").exists():
        st.warning("COMSOL kinematics files (`qdq_step.xlsx`, `dq.xlsx`) not found.")
    else:
        if st.button("▶️ Run Phase Lag Analysis", key="run_pl", type="primary"):
            with st.spinner("Computing forces with both velocity inputs..."):
                pl_results, pl_time = _run_phase_lag_analysis()
            st.session_state["pl_results"] = pl_results
            st.session_state["pl_time"] = pl_time

        if "pl_results" in st.session_state:
            pl_results = st.session_state["pl_results"]
            pl_time = st.session_state["pl_time"]

            # Metrics comparison table
            st.markdown("#### Metrics: Correct vs COMSOL-Matched Velocity")
            rows = []
            for phi in sorted(pl_results.keys()):
                for lbl in ["Heave", "Sway", "Surge"]:
                    mc = pl_results[phi]['metrics'][lbl]['correct']
                    mm = pl_results[phi]['metrics'][lbl]['match']
                    rows.append({
                        "Φ [°]": phi, "Component": lbl,
                        "Corr (correct)": f"{mc['corr']:.3f}",
                        "Lag (correct) [s]": f"{mc['lag']:+.3f}",
                        "NRMSE (correct) [%]": f"{mc['nrmse']:.1f}",
                        "Corr (matched)": f"{mm['corr']:.3f}",
                        "Lag (matched) [s]": f"{mm['lag']:+.3f}",
                        "NRMSE (matched) [%]": f"{mm['nrmse']:.1f}",
                    })
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

            # Phase lag detail plot for selected Φ
            phi_pl = st.selectbox("Select Φ for detailed view", [20, 30, 40, 50],
                                  index=3, key="pl_phi_select")
            r = pl_results[phi_pl]
            n = r['n']
            cfd = r['cfd']
            t = pl_time[:n]

            fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
            fig.suptitle(
                f"Phase Lag Analysis: Φ = {phi_pl}° — Effect of Joint Velocity Input\n"
                r"(Blue = correct $[\dot{q}_1,\dot{q}_2]$, Magenta = COMSOL-matched"
                r" $[\dot{q}_1,\dot{q}_1]$, Red = CFD)",
                fontsize=11)
            for k, (lbl, f_idx, c_idx) in enumerate(_FORCE_MAP):
                ax = axes[k]
                ax.plot(t, r['F_correct'][:, f_idx], 'b-', lw=1.5,
                        label=r'Analytical $[\dot{q}_1,\dot{q}_2]$')
                ax.plot(t, r['F_match'][:, f_idx], 'm:', lw=2,
                        label=r'Analytical $[\dot{q}_1,\dot{q}_1]$')
                ax.plot(t, cfd[:, c_idx], 'r--', lw=2, label='CFD (COMSOL)')
                ax.set_ylabel(f"{lbl} [N]")
                ax.legend(fontsize=9, loc='upper right')
                ax.grid(True, alpha=0.3)
                mc = r['metrics'][lbl]['correct']
                mm = r['metrics'][lbl]['match']
                ax.text(0.02, 0.92,
                        f"Correct: r={mc['corr']:.3f}, lag={mc['lag']:+.3f}s\n"
                        f"Matched: r={mm['corr']:.3f}, lag={mm['lag']:+.3f}s",
                        transform=ax.transAxes, fontsize=8, va='top',
                        bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', alpha=0.9))
            axes[-1].set_xlabel("Time [s]")
            fig.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

            # Summary bar charts
            phi_list = sorted(pl_results.keys())
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
            fig.suptitle(
                r"Impact of Joint Velocity: Correct $[\dot{q}_1,\dot{q}_2]$ vs"
                r" COMSOL-matched $[\dot{q}_1,\dot{q}_1]$", fontsize=11)
            x = np.arange(len(phi_list))
            w = 0.12
            for k, (lbl, _, _) in enumerate(_FORCE_MAP):
                corr_c = [pl_results[p]['metrics'][lbl]['correct']['corr'] for p in phi_list]
                corr_m = [pl_results[p]['metrics'][lbl]['match']['corr'] for p in phi_list]
                ax1.bar(x + k * w - w, corr_c, w, label=f'{lbl} (correct)', alpha=0.6)
                ax1.bar(x + k * w - w + 3 * w, corr_m, w, label=f'{lbl} (matched)', alpha=0.9)
            ax1.set_xticks(x + w)
            ax1.set_xticklabels([f"{p}°" for p in phi_list])
            ax1.set_ylabel("Correlation with CFD")
            ax1.set_ylim(0.6, 1.05)
            ax1.legend(fontsize=7, ncol=2)
            ax1.grid(True, alpha=0.3, axis='y')
            ax1.set_title("Pearson Correlation")

            for k, (lbl, _, _) in enumerate(_FORCE_MAP):
                lag_c = [pl_results[p]['metrics'][lbl]['correct']['lag'] for p in phi_list]
                lag_m = [pl_results[p]['metrics'][lbl]['match']['lag'] for p in phi_list]
                ax2.bar(x + k * w - w, lag_c, w, label=f'{lbl} (correct)', alpha=0.6)
                ax2.bar(x + k * w - w + 3 * w, lag_m, w, label=f'{lbl} (matched)', alpha=0.9)
            ax2.set_xticks(x + w)
            ax2.set_xticklabels([f"{p}°" for p in phi_list])
            ax2.set_ylabel("Phase Lag [s]")
            ax2.legend(fontsize=7, ncol=2)
            ax2.grid(True, alpha=0.3, axis='y')
            ax2.axhline(y=0, color='k', lw=0.5)
            ax2.set_title("Phase Lag (negative = CFD lags)")
            fig.tight_layout()
            st.pyplot(fig)
            plt.close(fig)

            # Key Finding
            st.divider()
            st.subheader("Key Findings")

            avg_corr_c = np.mean([pl_results[p]['metrics']['Heave']['correct']['corr'] for p in phi_list])
            avg_corr_m = np.mean([pl_results[p]['metrics']['Heave']['match']['corr'] for p in phi_list])
            avg_lag_c = np.mean([abs(pl_results[p]['metrics']['Heave']['correct']['lag']) for p in phi_list])
            avg_lag_m = np.mean([abs(pl_results[p]['metrics']['Heave']['match']['lag']) for p in phi_list])

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Heave Corr. (correct)", f"{avg_corr_c:.3f}")
            c2.metric("Heave Corr. (matched)", f"{avg_corr_m:.3f}", f"{avg_corr_m - avg_corr_c:+.3f}")
            c3.metric("Avg |Phase Lag| (correct)", f"{avg_lag_c:.3f} s")
            c4.metric("Avg |Phase Lag| (matched)", f"{avg_lag_m:.3f} s", f"{avg_lag_m - avg_lag_c:+.3f} s", delta_color="inverse")

            st.markdown(f"""
#### Root Cause of CFD–Analytical Discrepancy

1. **COMSOL's moving mesh** derived membrane velocity only from **joint 1** ($\\dot{{q}}_1$),
   effectively applying $[\\dot{{q}}_1, \\dot{{q}}_1]$ instead of $[\\dot{{q}}_1, \\dot{{q}}_2]$.

2. **When the analytical model uses the same velocity limitation**, correlation with CFD
   improves from **{avg_corr_c:.3f} → {avg_corr_m:.3f}** (Heave average) and the phase lag
   drops from **{avg_lag_c:.3f}s → {avg_lag_m:.3f}s**.

3. The **remaining ~20–30% NRMSE** is attributable to the fundamental difference between
   quasi-steady Resistive Force Theory (RFT) and the full Navier–Stokes solution in COMSOL.

4. The `[dq1, dq1]` pattern found in the original MATLAB comparison script was **not a bug** —
   it was intentionally matching the COMSOL model's velocity limitation for a fair comparison.
""")


def tab_matlab_corrections():
    """MATLAB Bug Documentation — corrections identified during reimplementation."""
    st.header("Corrections Applied to the Original MATLAB Code")
    st.caption("Bugs identified and fixed during the Python reimplementation of the MATLAB codebase.")

    st.markdown("""
During the Python reimplementation, several errors were identified in the original MATLAB
codebase (`undFinDynamics.m`, `hydrodynamicforces_Func.m`, `hydrodynamictorques_Func.m`).
These bugs were verified by symbolic comparison and point-by-point numerical testing against
the original MATLAB Symbolic Math Toolbox output. Below is a summary of each issue and
how it was resolved.
""")

    # Bug 1: vn sign
    with st.expander("**Bug 1 — Normal Velocity Sign Error (vn quadratic term)**", expanded=True):
        st.markdown(r"""
**Location:** `hydrodynamicforces_Func.m` / `hydrodynamictorques_Func.m` — normal velocity $v_n$ computation (Eq. 14)

**Description:**
The normal velocity formula contains a quadratic term in $\bar{w}$ (the normalized membrane coordinate).
The MATLAB Symbolic Toolbox generates the expression:

$$v_n = \frac{Dh}{\delta}\left[\dot{q}_1 + (1-\cos\Delta q)(\dot{q}_1+\dot{q}_2)\bar{w}^2
       + (\cos\Delta q\,(\dot{q}_1+\dot{q}_2) - 2\dot{q}_1)\bar{w}\right]$$

In the original Python port (prior to correction), the $\bar{w}^2$ term was coded with a
**negative** sign:

```python
# INCORRECT (original):
vn_inner = dq1 - (1-cos_dq)*(dq1+dq2)*w_bar**2 + ...
```

The MATLAB symbolic output, after simplifying the denominator
$h^2(\cos\Delta q - 1)\cdot 2 - D^2 = -w_\mathrm{max}^2$, yields a **positive** sign:

```python
# CORRECT (fixed):
vn_inner = dq1 + (1-cos_dq)*(dq1+dq2)*w_bar**2 + ...
```

**Impact:** This sign error propagates into the $v_n^2$ drag formulation, affecting all three
force components (Heave, Sway, Surge) and the hydrodynamic torques. The overall force magnitude
and phase shift are distorted at every time step.

**Verification:** After fixing, the Python integrand matches the MATLAB symbolic expression
to <$10^{-12}$ relative error across 375 test points (5 operating points × 5 heights × 5 widths
× 3 components).
""")

    # Bug 2: Elastic torques
    with st.expander("**Bug 2 — Elastic Torque Panel Assignment (N > 2 rays)**", expanded=True):
        st.markdown(r"""
**Location:** `undFinDynamics.m` — elastic torque loop

**Description:**
The MATLAB code computes elastic membrane torques $\tau_G$ with the following loop structure:

```matlab
for i = 1:N
    if i < N
        j = 1;       % panel 1-2
    else
        j = 2;       % panel 2-3
    end
    tauG(:,i) = (-1)^(i+1) * sin(q(j)-q(j+1)) * EK * ...
                 (DeltaFunc(D,hmax,q([j:j+1])) - DeltaFunc(D,hmin,q([j:j+1])));
end
```

For **N = 3 rays**, the loop assigns:
| Ray (i) | j | Panel Used | Correct Panel |
|:-------:|:-:|:----------:|:------------:|
| 1 | 1 | Panel 1→2 | Panel 1→2 ✅ |
| 2 | 1 | Panel 1→2 | Panel **1→2 and 2→3** ❌ |
| 3 | 2 | Panel 2→3 | Panel 2→3 ✅ |

**Problem:** Ray 2 (interior ray) receives its elastic torque only from panel 1→2, but
the correct physics requires contributions from **both** adjacent membranes (panel 1→2
and panel 2→3). This logic works correctly only for N = 2 rays.

**Impact:** For N ≥ 3, interior rays experience incorrect elastic restoring forces.
At N = 2, the bug is dormant (no interior rays).

**Fix:** The Python implementation sums contributions from both adjacent panels for each
interior ray.
""")

    # Bug 3: Hydrodynamic torque missing contribution
    with st.expander("**Bug 3 — Missing Hydrodynamic Torque Contribution (N > 2)**", expanded=True):
        st.markdown(r"""
**Location:** `hydrodynamictorques_Func.m` — torque return vector

**Description:**
The MATLAB torque function returns a vector $T = [T_1, T_2, T_3]$ for N = 3 rays:
- $T_1$ = torque on joint 1 from panel 1→2
- $T_2$ = torque on joint 2 from panel 1→2 (far-end lever arm $w/w_\mathrm{max}$)
- $T_3$ = torque on joint 3 from panel 2→3 (far-end lever arm $w/w_\mathrm{max}$)

**Problem:** Joint 2 also receives a torque contribution from panel 2→3 (where it acts as
the near-end joint), but this contribution is **not included** in the return vector. The
MATLAB function only returns one torque per joint, even though interior joints border two
panels.

**Impact:** For N ≥ 3, interior rays are missing part of their hydrodynamic torque input.
Similar to Bug 2, this is invisible at N = 2.

**Fix:** The Python implementation accumulates both near-end and far-end torque contributions
for every joint from all adjacent panels.
""")

    # Bug 4: Force component labels
    with st.expander("**Bug 4 — Force Component Label Swapping**", expanded=True):
        st.markdown(r"""
**Location:** `undFinDynamics.m` — plotting section

**Description:**
The MATLAB integrand computes three force components from the membrane surface integral:
- $F_x$: uses $\sin(q)$ terms (lateral direction)
- $F_y$: uses $\cos(q)$ terms (vertical direction)
- $F_z$: uses $\sin(\Delta q) \cdot h/D$ terms (propulsive direction)

The plotting code in `undFinDynamics.m` uses the legend:
```matlab
legend('F_{heave}', 'F_{sway}', 'F_{surge}')
```
mapping $F(1) = F_x \to$ "Heave" and $F(2) = F_y \to$ "Sway".

**Problem:** Physically, Heave (vertical/lift) corresponds to $\cos(q)$ terms = $F_y$,
and Sway (lateral) corresponds to $\sin(q)$ terms = $F_x$. The MATLAB code swaps these
labels.

**Correct mapping:**

| Python Index | Math Terms | Physical Meaning | MATLAB Label (incorrect) |
|:---:|:---:|:---:|:---:|
| F[0] | sin(q) terms | **Sway** (lateral) | "Heave" |
| F[1] | cos(q) terms | **Heave** (vertical) | "Sway" |
| F[2] | sin(Δq) terms | **Surge** (propulsive) | "Surge" ✅ |

**Impact:** Does not affect numerical results, but produces misleading plots and
makes CFD validation confusing. The Python implementation uses the physically correct
labels and applies the proper column mapping when comparing with COMSOL data.
""")

    st.divider()
    st.markdown("""
> **Note:** Bugs 2 and 3 only manifest for **N ≥ 3 rays**. Since the COMSOL CFD
> validation was done with **N = 2 rays**, these bugs did not affect the published
> comparison results. Bug 1 (vn sign) and Bug 4 (label swap) affect all
> configurations.
""")


def tab_about():
    """About / reference page."""
    st.header("About This Application")
    st.markdown("""
### Kinematic and Dynamic Modeling of Biomimetic Undulating Fins

This application implements the analytical model from the Master's thesis:

> **Alhamzah Ihsan Ali**, *"Kinematic and Dynamic Modeling of Biomimetic Undulating Fins 
> for Underwater Robots"*, Bahçeşehir University, Jan 2021.  
> Supervisor: Assoc. Prof. Dr. Berke Gür

---

### 📖 Equations Reference

| Equation | Description | Tab |
|----------|-------------|-----|
| Eq. (1)–(3) | Traveling wave kinematics | Kinematics |
| Eq. (4)–(7) | Lagrangian dynamics: τ = L + Γ + Φ | Time Domain |
| Eq. (8)–(10) | Elastic membrane deflection & torques | Time Domain |
| Eq. (12)–(16) | Hydrodynamic forces via resistive theory | Time Domain |
| Ch. 5 | Parametric studies (A, f, Φ, D, N) | Parametric |
| Eq. (18)–(19) | Dual-fin robot body dynamics | Robot |

---

### 🔬 Model Overview

The undulating fin consists of **N** rigid rays connected by elastic membranes, 
driven by sinusoidal oscillation with inter-ray phase shift **Φ** to create a 
traveling wave along the fin. The model computes:

1. **Kinematics**: Prescribed ray oscillation angles
2. **Hydrodynamics**: Membrane drag forces via quasi-steady resistive force theory
3. **Dynamics**: Full Lagrangian formulation including inertia, friction, gravity, 
   membrane elasticity, and hydrodynamic loads
4. **Robot motion**: Dual-fin underwater robot trajectory integration

---

### 🛠️ Technology Stack

- **Python 3.10+** with NumPy, SciPy, Matplotlib
- **Streamlit** for the interactive dashboard
- **scipy.integrate.dblquad** for double integration over membrane surfaces
""")



# ──────────────────────────── OpenFOAM Pipeline ─────────────────────

_OF_BASE = Path(__file__).resolve().parent / "openfoam_cases"


def tab_openfoam(params: FinParams):
    """OpenFOAM validation pipeline tab."""
    st.header("🔬 OpenFOAM Validation Pipeline")
    st.caption(
        "Generate, manage, and validate OpenFOAM FSI simulations "
        "against analytical models (RFT / Morison / Lighthill)."
    )

    sec1, sec2, sec3 = st.tabs([
        "⚙️ Case Generator",
        "📂 Results Browser",
        "📊 Validation",
    ])

    # ── Section 1: Case Generator ──────────────────────────────────
    with sec1:
        _of_case_generator(params)

    # ── Section 2: Results Browser ─────────────────────────────────
    with sec2:
        _of_results_browser()

    # ── Section 3: Validation ──────────────────────────────────────
    with sec3:
        _of_validation(params)


def _of_case_generator(params: FinParams):
    st.subheader("Generate OpenFOAM Case")
    st.markdown("""
    Creates a ready-to-run case directory with:
    - **blockMeshDict** — parametric channel mesh
    - **Boundary conditions** — movingWallVelocity on fin patch
    - **dynamicMeshDict** — ALE morphing via `displacementLaplacian`
    - **controlDict** — force function objects for data extraction
    - **fin_motion.py** — kinematics script driven by FinParams
    """)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Simulation Settings**")
        n_cycles = st.number_input("Number of cycles", 1, 20, 3, key="of_ncyc")
        writes_per_cycle = st.number_input("Write intervals / cycle", 5, 100, 20, key="of_wpc")
        max_co = st.number_input("Max Courant number", 0.1, 5.0, 1.0, 0.1, key="of_maxco")
        solver = st.selectbox("Solver", ["pimpleDyMFoam", "overPimpleDyMFoam", "pimpleFoam"],
                              key="of_solver")

    with col2:
        st.markdown("**Domain & Mesh**")
        cells_per_D = st.slider("Cells per D", 10, 100, 40, key="of_cpd")
        inlet_len = st.number_input("Inlet length [×D]", 1.0, 20.0, 5.0, 1.0, key="of_inlet")
        outlet_len = st.number_input("Outlet length [×D]", 5.0, 30.0, 15.0, 1.0, key="of_outlet")
        lateral = st.number_input("Lateral width [×D]", 2.0, 20.0, 5.0, 1.0, key="of_lat")
        nu = st.number_input("Kinematic viscosity ν [m²/s]", 1e-7, 1e-3, 1e-6, format="%e",
                             key="of_nu")

    case_name = st.text_input(
        "Case directory name",
        value=f"phi{int(np.rad2deg(params.Phi))}_A{int(np.rad2deg(params.A))}"
              f"_f{params.f:.0f}_N{params.N}",
        key="of_casename",
    )

    # Show fin parameters summary
    with st.expander("📋 Fin parameters for this case"):
        pcols = st.columns(4)
        pcols[0].metric("N (rays)", params.N)
        pcols[1].metric("A [°]", f"{np.rad2deg(params.A):.1f}")
        pcols[2].metric("f [Hz]", f"{params.f:.1f}")
        pcols[3].metric("Φ [°]", f"{np.rad2deg(params.Phi):.1f}")
        pcols2 = st.columns(4)
        pcols2[0].metric("D [m]", f"{params.D:.3f}")
        pcols2[1].metric("hmax [m]", f"{params.hmax:.3f}")
        pcols2[2].metric("hmin [m]", f"{params.hmin:.3f}")
        pcols2[3].metric("ρ [kg/m³]", f"{params.rho:.0f}")

    case_path = _OF_BASE / case_name
    overwrite = st.checkbox("Overwrite if exists", value=False, key="of_overwrite")

    if st.button("🚀 Generate Case", type="primary", key="of_gen_btn"):
        cfg = OFConfig(
            n_cycles=n_cycles,
            write_interval_per_cycle=writes_per_cycle,
            max_courant=max_co,
            solver=solver,
            cells_per_D=cells_per_D,
            inlet_length=inlet_len,
            outlet_length=outlet_len,
            lateral_width=lateral,
            nu=nu,
        )
        try:
            result_path = generate_case(params, cfg, case_path, overwrite=overwrite)
            st.success(f"✅ Case generated at: `{result_path}`")
            st.info("To run: `cd " + str(result_path) + " && ./Allrun`")

            # Show directory structure
            with st.expander("📁 Generated files"):
                _show_tree(result_path)
        except FileExistsError:
            st.error(
                f"Case directory already exists: `{case_path}`. "
                "Enable *Overwrite if exists* or choose a different name."
            )
        except Exception as e:
            st.error(f"Error generating case: {e}")

    # Batch generation
    st.divider()
    st.subheader("Batch Case Generation")
    st.markdown("Generate cases for multiple phase shifts at once (for parametric study).")

    phi_values_str = st.text_input(
        "Phase shifts Φ [°] (comma-separated)",
        value="0, 20, 30, 40, 50",
        key="of_batch_phi",
    )

    if st.button("🔄 Generate Batch", key="of_batch_btn"):
        try:
            phi_list = [float(x.strip()) for x in phi_values_str.split(",")]
        except ValueError:
            st.error("Invalid input — enter comma-separated numbers.")
            return

        cfg = OFConfig(
            n_cycles=n_cycles,
            write_interval_per_cycle=writes_per_cycle,
            max_courant=max_co,
            solver=solver,
            cells_per_D=cells_per_D,
            inlet_length=inlet_len,
            outlet_length=outlet_len,
            lateral_width=lateral,
            nu=nu,
        )
        progress = st.progress(0, text="Generating cases...")
        generated = []
        for i, phi_deg in enumerate(phi_list):
            p = FinParams(
                N=params.N, A=params.A, f=params.f,
                Phi=np.deg2rad(phi_deg), X=params.X, D=params.D,
                hmax=params.hmax, hmin=params.hmin,
                Cn=params.Cn, Cm=params.Cm, rho=params.rho,
            )
            cname = f"phi{int(phi_deg)}_A{int(np.rad2deg(params.A))}_f{params.f:.0f}_N{params.N}"
            cpath = _OF_BASE / cname
            try:
                generate_case(p, cfg, cpath, overwrite=overwrite)
                generated.append(cname)
            except FileExistsError:
                st.warning(f"Skipped (exists): {cname}")
            progress.progress((i + 1) / len(phi_list), text=f"Generated {i+1}/{len(phi_list)}")

        if generated:
            st.success(f"✅ Generated {len(generated)} cases: {', '.join(generated)}")


def _show_tree(path: Path, prefix: str = "", max_depth: int = 3, _depth: int = 0):
    """Display directory tree in st.code."""
    if _depth >= max_depth:
        return
    entries = sorted(path.iterdir(), key=lambda p: (not p.is_dir(), p.name))
    lines = []
    for entry in entries:
        if entry.is_dir():
            lines.append(f"{prefix}📁 {entry.name}/")
            sub = sorted(entry.iterdir(), key=lambda p: (not p.is_dir(), p.name))
            for s in sub:
                lines.append(f"{prefix}    {'📁' if s.is_dir() else '📄'} {s.name}")
        else:
            lines.append(f"{prefix}📄 {entry.name}")
    st.code("\n".join(lines), language=None)


def _of_results_browser():
    st.subheader("Completed Simulations")

    if not _OF_BASE.is_dir():
        st.info("No OpenFOAM cases found yet. Generate cases in the *Case Generator* tab.")
        return

    # Also allow loading from a custom directory
    custom_dir = st.text_input(
        "Cases directory",
        value=str(_OF_BASE),
        key="of_results_dir",
    )
    search_dir = Path(custom_dir)

    if not search_dir.is_dir():
        st.warning(f"Directory not found: {search_dir}")
        return

    # List all subdirectories
    all_dirs = sorted([d for d in search_dir.iterdir() if d.is_dir()])
    if not all_dirs:
        st.info("No case directories found.")
        return

    st.markdown(f"Found **{len(all_dirs)}** case directories:")

    rows = []
    for d in all_dirs:
        has_mesh = (d / "constant" / "polyMesh").is_dir() or (d / "system" / "blockMeshDict").is_file()
        has_results = False
        try:
            find_force_file(d)
            has_results = True
        except FileNotFoundError:
            pass

        status = "✅ Results" if has_results else ("🔧 Meshed" if has_mesh else "📝 Generated")
        rows.append({
            "Case": d.name,
            "Status": status,
            "Path": str(d),
        })

    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    # Preview a specific case
    completed = list_completed_cases(search_dir)
    if completed:
        st.divider()
        st.subheader("Preview Force Data")
        sel = st.selectbox(
            "Select completed case",
            [c.name for c in completed],
            key="of_preview_case",
        )
        sel_path = search_dir / sel
        if st.button("📈 Load Forces", key="of_load_btn"):
            try:
                fd = load_forces(sel_path)
                st.success(f"Loaded {len(fd.time)} time steps ({fd.time[0]:.4f} – {fd.time[-1]:.4f} s)")

                fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
                labels = ["Fx (Heave)", "Fy (Sway)", "Fz (Surge)"]
                for k in range(3):
                    axes[k].plot(fd.time, fd.total[:, k], "b-", lw=1, label="Total")
                    axes[k].plot(fd.time, fd.pressure[:, k], "r--", lw=0.8, label="Pressure")
                    axes[k].plot(fd.time, fd.viscous[:, k], "g:", lw=0.8, label="Viscous")
                    axes[k].set_ylabel(f"{labels[k]} [N]")
                    axes[k].legend(fontsize=8)
                    axes[k].grid(True, alpha=0.3)
                axes[-1].set_xlabel("Time [s]")
                fig.suptitle(f"OpenFOAM Forces — {sel}", fontsize=13)
                fig.tight_layout()
                st.pyplot(fig)
                plt.close(fig)
            except Exception as e:
                st.error(f"Error loading forces: {e}")
    else:
        st.info("No completed simulations found (no postProcessing data).")


def _of_validation(params: FinParams):
    st.subheader("OpenFOAM vs. Analytical Validation")

    st.markdown("""
    Compare OpenFOAM CFD results against the analytical models.
    Uses the same metrics as the COMSOL comparison:
    **RMSE**, **NRMSE**, and **Pearson correlation**.
    """)

    # Let user pick a case directory
    custom_dir = st.text_input(
        "Cases directory",
        value=str(_OF_BASE),
        key="of_val_dir",
    )
    search_dir = Path(custom_dir)
    completed = list_completed_cases(search_dir) if search_dir.is_dir() else []

    if not completed:
        st.info(
            "No completed OpenFOAM cases found. Run simulations first, "
            "or point to a directory containing completed cases."
        )
        return

    sel_names = st.multiselect(
        "Select cases to validate",
        [c.name for c in completed],
        default=[c.name for c in completed[:3]],
        key="of_val_cases",
    )

    methods = st.multiselect(
        "Analytical methods",
        ["rft", "morison", "lighthill"],
        default=["rft", "morison"],
        key="of_val_methods",
    )

    # Component mapping
    with st.expander("🔧 Component Mapping"):
        st.markdown("""
        Map analytical force components to OpenFOAM output columns.
        Default: `Heave→Fx, Sway→Fy, Surge→Fz`.
        """)
        cx = st.selectbox("Heave → OF column", [0, 1, 2], index=0, key="of_cx")
        cy = st.selectbox("Sway → OF column", [0, 1, 2], index=1, key="of_cy")
        cz = st.selectbox("Surge → OF column", [0, 1, 2], index=2, key="of_cz")
        cmap = {0: cx, 1: cy, 2: cz}

    phi_override = st.number_input(
        "Override Φ [°] (used for analytical model)",
        value=float(np.rad2deg(params.Phi)),
        key="of_phi",
    )

    if st.button("🔍 Run Validation", type="primary", key="of_val_btn"):
        val_params = FinParams(
            N=params.N, A=params.A, f=params.f,
            Phi=np.deg2rad(phi_override), X=params.X, D=params.D,
            hmax=params.hmax, hmin=params.hmin,
            Cn=params.Cn, Cm=params.Cm, rho=params.rho,
        )

        all_results = []
        progress = st.progress(0, text="Validating...")

        for i, name in enumerate(sel_names):
            case_path = search_dir / name
            try:
                results = validate_case(
                    case_path, val_params,
                    methods=methods, component_map=cmap,
                )
                all_results.extend(results)
            except Exception as e:
                st.warning(f"Error validating {name}: {e}")
            progress.progress((i + 1) / len(sel_names))

        if not all_results:
            st.error("No validation results produced.")
            return

        # Metrics table
        st.divider()
        st.subheader("📊 Metrics Summary")
        df = results_to_dataframe(all_results)
        st.dataframe(df, use_container_width=True, hide_index=True)

        # Best method per case
        best = df.groupby(["Case", "Method"])["NRMSE [%]"].mean().reset_index()
        best_idx = best.groupby("Case")["NRMSE [%]"].idxmin()
        st.markdown("**Best method per case (lowest mean NRMSE):**")
        st.dataframe(best.loc[best_idx], use_container_width=True, hide_index=True)

        # Force overlay plots
        st.divider()
        st.subheader("📈 Force Comparison Plots")

        for r in all_results:
            with st.expander(f"{r.case_name} — {r.method.upper()} (NRMSE: {r.mean_nrmse:.1f}%)"):
                fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
                for k, label in enumerate(FORCE_LABELS):
                    axes[k].plot(r.time, r.F_of[:, k], "b-", lw=1.2, label="OpenFOAM")
                    axes[k].plot(r.time, r.F_an[:, k], "r--", lw=1.2, label=r.method.upper())
                    axes[k].set_ylabel(f"{label} [N]")
                    axes[k].legend(fontsize=8)
                    axes[k].grid(True, alpha=0.3)
                    m = r.metrics[k]
                    axes[k].set_title(
                        f"{label}: NRMSE={m.nrmse_pct:.1f}%, "
                        f"Corr={m.correlation:.3f}", fontsize=10,
                    )
                axes[-1].set_xlabel("Time [s]")
                fig.suptitle(f"{r.case_name} — {r.method.upper()}", fontsize=13)
                fig.tight_layout()
                st.pyplot(fig)
                plt.close(fig)

        # NRMSE bar chart
        st.divider()
        st.subheader("📊 NRMSE Comparison")
        fig, ax = plt.subplots(figsize=(10, 4))
        method_names = sorted(set(r.method for r in all_results))
        x = np.arange(len(FORCE_LABELS))
        width = 0.8 / max(len(method_names), 1)
        for j, meth in enumerate(method_names):
            matching = [r for r in all_results if r.method == meth]
            if matching:
                nrmses = [matching[0].metrics[k].nrmse_pct for k in range(3)]
                ax.bar(x + j * width, nrmses, width, label=meth.upper(), alpha=0.8)
        ax.set_xticks(x + width * (len(method_names) - 1) / 2)
        ax.set_xticklabels(FORCE_LABELS)
        ax.set_ylabel("NRMSE [%]")
        ax.set_title("Normalized RMSE by Component")
        ax.legend()
        ax.grid(True, alpha=0.3, axis="y")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)


# ──────────────────────────── Mathematical Models ────────────────────


def tab_math_models():
    """Comprehensive documentation of all mathematical models."""
    st.header("📐 Mathematical Models Reference")
    st.caption(
        "Complete equation reference for all models implemented in this application. "
        "Organized by module for quick lookup."
    )

    # ─────────── 1. Kinematics ───────────
    st.subheader("1. Traveling-Wave Kinematics")
    st.markdown("""
Each of the $N$ fin rays oscillates sinusoidally with an inter-ray phase shift
$\\Phi$ that generates a traveling wave along the fin (Thesis Sec. 2.2).
""")

    st.markdown("**Eq. (1)** — Joint angle:")
    st.latex(r"""
q_i(t) = A \sin\!\bigl(2\pi f t - (i-1)\Phi\bigr) + X
\qquad i = 1, \dots, N
""")

    st.markdown("**Eq. (2)** — Angular velocity:")
    st.latex(r"""
\dot{q}_i(t) = 2\pi f A \cos\!\bigl(2\pi f t - (i-1)\Phi\bigr)
""")

    st.markdown("**Eq. (2b)** — Angular acceleration:")
    st.latex(r"""
\ddot{q}_i(t) = -\bigl(2\pi f\bigr)^2 A \sin\!\bigl(2\pi f t - (i-1)\Phi\bigr)
""")

    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $A$ | Oscillation amplitude | rad |
| $f$ | Frequency | Hz |
| $\\Phi$ | Inter-ray phase shift | rad |
| $X$ | Bias angle | rad |
| $N$ | Number of rays | — |
""")

    st.markdown("**Eq. (3)** — Derived quantities:")
    st.latex(r"""
\lambda = \frac{2\pi D}{|\Phi|}
\qquad
w_n = \frac{|\Phi|(N-1)}{2\pi}
""")
    st.caption("$\\lambda$: wavelength [m], $w_n$: number of waves on fin.")

    st.divider()

    # ─────────── 2. Dynamics ───────────
    st.subheader("2. Lagrangian Dynamics")
    st.markdown("**Eq. (7)** — Equation of motion for each ray joint:")
    st.latex(r"""
\tau_i = \underbrace{M_i \ddot{q}_i + b_i \dot{q}_i + k_{c,i}
\operatorname{sgn}(\dot{q}_i) + g_i \sin(q_i)}_{L_i}
+ \underbrace{\Gamma_i(\Delta q)}_{\text{Elastic}}
+ \underbrace{\Phi_i(q, \dot{q})}_{\text{Hydro}}
""")

    st.markdown("where each term represents:")
    st.markdown("""
| Term | Expression | Description |
|------|------------|-------------|
| Inertia | $M_i \\ddot{q}_i$ | $M_i = m_i (l_i/2)^2 + I_{z,i}$ |
| Viscous friction | $b_i \\dot{q}_i$ | Linear damping |
| Coulomb friction | $k_{c,i} \\operatorname{sgn}(\\dot{q}_i)$ | Dry friction |
| Gravity | $g_i \\sin(q_i)$ | $g_i = (l_i/2) m_i g$ |
| Elastic (membrane) | $\\Gamma_i(\\Delta q)$ | See §3 below |
| Hydrodynamic | $\\Phi_i(q, \\dot{q})$ | See §4–6 below |
""")

    st.divider()

    # ─────────── 3. Elastic Membrane ───────────
    st.subheader("3. Elastic Membrane Torques")
    st.markdown("""
The elastic membrane connecting rays $j$ and $j+1$ exerts restoring torques
(Thesis Eq. 8–10). The deflection function $\\Lambda(x)$ captures geometric
nonlinearity:
""")

    st.markdown("**Eq. (10)** — Deflection function:")
    st.latex(r"""
\Lambda(x) = \frac{x^3 \alpha}{3D}
\left(1 + \frac{\nu\bigl(8 - 3\alpha\sqrt{\kappa^2+1}\bigr)}{4}\right)
- \frac{D^2}{2p^3}
\left(\sqrt{\kappa^4+\kappa^2} - \operatorname{sinh}^{-1}(\kappa)\right)
\left(1 + \frac{\nu(4+\alpha^2)}{4}\right)
""")

    st.markdown("where:")
    st.latex(r"""
\alpha = \frac{D}{\ell}, \qquad
p = \sqrt{2(1-\cos\Delta q)}, \qquad
\kappa = \frac{x}{D}\,p
""")

    st.markdown("""
The torque contribution from membrane $\\Omega_{j|j+1}$ on joint $j$ is:
""")
    st.markdown("**Eq. (8–9)** — Membrane torque on joint $j$:")
    st.latex(r"""
\tau_{\Gamma,j} = (-1)^{j+1} \sin(q_j - q_{j+1})\; E_K \;
\bigl[\Lambda(h_{\max}) - \Lambda(h_{\min})\bigr]
""")

    st.divider()

    # ─────────── 4. RFT ───────────
    st.subheader("4. Model I — Resistive Force Theory (RFT)")
    st.markdown("""
The original quasi-steady drag model from the thesis. Forces are purely
resistive — proportional to $v_n^2$ (Thesis Sec. 2.3.2).
""")

    st.markdown("**Eq. (12)** — Force per unit area on the membrane surface:")
    st.latex(r"""
\mathbf{f}_n = -\tfrac{1}{2}\,\rho\, C_n\, \mathbf{e}_n\, v_n^2\,
\operatorname{sgn}(v_n)
""")

    st.markdown("**Eq. (13)** — Total force — double integration over each panel $\\Omega_{i|i+1}$:")
    st.latex(r"""
\mathbf{F}_N = \sum_{i=1}^{N-1}
\int_{h_{\min}}^{h_{\max}} \int_0^{w_{\max}(h)}
\mathbf{f}_n \; dw\; dh
""")

    st.markdown("**Eq. (14)** — Normal velocity at membrane point $(h, \\bar{w})$:")
    st.latex(r"""
v_n = \frac{Dh}{\delta}\Bigl[
\dot{q}_1 + (1-\cos\Delta q)(\dot{q}_1+\dot{q}_2)\bar{w}^2
+ \bigl(\cos\Delta q\,(\dot{q}_1+\dot{q}_2) - 2\dot{q}_1\bigr)\bar{w}
\Bigr]
""")

    st.markdown("**Eq. (15)** — Geometric factor:")
    st.latex(r"""
\delta = \sqrt{2D^2\bar{w}(\bar{w}-1)(1-\cos\Delta q) + h^2\sin^2\!\Delta q + D^2}
""")

    st.markdown("**Eq. (16)** — Normal vector:")
    st.latex(r"""
\mathbf{e}_n = \frac{D}{\delta}
\begin{bmatrix}
\sin q_1\,(\bar{w}-1) - \sin q_2\,\bar{w} \\
-\cos q_1\,(\bar{w}-1) + \cos q_2\,\bar{w} \\
-\sin\Delta q\; h/D
\end{bmatrix}
""")

    st.markdown("**Integration limit** (strip width at height $h$):")
    st.latex(r"""
w_{\max}(h) = \sqrt{2h^2(1-\cos\Delta q) + D^2}
""")

    with st.expander("Parameter summary"):
        st.markdown("""
| Symbol | Description | Typical value |
|--------|-------------|---------------|
| $C_n$ | Normal drag coefficient | 2.65 – 2.7 |
| $\\rho$ | Fluid density | 1000 – 1022 kg/m³ |
| $D$ | Inter-ray distance | 0.09 – 0.10 m |
| $h_{\\min}, h_{\\max}$ | Membrane height bounds | 0 – 0.175 m |
| $\\Delta q$ | Angle difference $q_i - q_{i+1}$ | — |
| $\\bar{w}$ | Normalized coordinate $w/w_{\\max}$ | [0, 1] |
""")

    st.divider()

    # ─────────── 5. Morison ───────────
    st.subheader("5. Model II — Morison's Equation (Drag + Inertia)")
    st.markdown("""
Extends RFT by adding an **added-mass inertia term** (Morison et al., 1950).
This accounts for the reaction force of the fluid accelerated by the moving membrane.
""")

    st.markdown("**Eq. (M1)** — Force per unit area:")
    st.latex(r"""
\mathbf{f}_n = \underbrace{-\tfrac{1}{2}\,\rho\, C_n\, \mathbf{e}_n\, v_n^2\,
\operatorname{sgn}(v_n)}_{\text{Drag (RFT)}}
\;\underbrace{-\;\rho\, C_m\, \frac{\pi}{4}\, w_{\max}(h)\;
\mathbf{e}_n\, \dot{v}_n}_{\text{Added-mass inertia}}
""")

    st.markdown("""
**Physical interpretation:**
- The drag term resists motion — proportional to $v_n^2$
- The inertia term resists **acceleration** — proportional to $\\dot{v}_n$
- Both forces are **negative** (reaction forces on the body moving in still fluid)
""")

    st.markdown("**Added mass model** — for a thin plate of chord $c$ oscillating normally:")
    st.latex(r"""
m_a = \rho\,\frac{\pi}{4}\, c^2
\quad \text{(per unit span)}
\qquad \Rightarrow \qquad
\frac{m_a}{c} = \rho\,\frac{\pi}{4}\, c
\quad \text{(per unit area)}
""")
    st.caption("Here $c = w_{\\max}(h)$ is the local membrane chord at height $h$.")

    st.markdown("**Eq. (M2)** — Normal acceleration $\\dot{v}_n$ — analytic time derivative of Eq. (14):")
    st.latex(r"""
\dot{v}_n = \frac{Dh}{\delta}\Bigl[
\dot{A} + \dot{B}\,\bar{w}^2 + \dot{C}\,\bar{w}
\Bigr]
- \frac{Dh\,\dot{\delta}}{\delta^2}\Bigl[
A + B\,\bar{w}^2 + C\,\bar{w}
\Bigr]
""")

    with st.expander("Coefficient derivatives (chain rule)"):
        st.latex(r"""
\begin{aligned}
A &= \dot{q}_1, &
\dot{A} &= \ddot{q}_1 \\[4pt]
B &= (1-\cos\Delta q)(\dot{q}_1+\dot{q}_2), &
\dot{B} &= \sin\Delta q\;\delta_q\;(\dot{q}_1+\dot{q}_2)
+ (1-\cos\Delta q)(\ddot{q}_1+\ddot{q}_2) \\[4pt]
C &= \cos\Delta q\;(\dot{q}_1+\dot{q}_2) - 2\dot{q}_1, &
\dot{C} &= -\sin\Delta q\;\delta_q\;(\dot{q}_1+\dot{q}_2)
+ \cos\Delta q\;(\ddot{q}_1+\ddot{q}_2) - 2\ddot{q}_1
\end{aligned}
""")
        st.markdown("where $\\delta_q = \\dot{q}_1 - \\dot{q}_2$ is the angular velocity difference.")
        st.latex(r"""
\dot{\delta}^2 = 2D^2\bar{w}(\bar{w}-1)\sin\Delta q\;\delta_q
+ 2h^2\sin\Delta q\cos\Delta q\;\delta_q
""")

    with st.expander("Cm tuning guidance"):
        st.markdown("""
| $C_m$ value | Meaning |
|-------------|---------|
| 1.0 | Theoretical flat plate (ideal) |
| 0.5 – 0.8 | Typical for flexible membranes |
| 0.0 | Reduces to pure RFT (no inertia) |

**How to tune:** Use the **CFD Comparison → Morison Cm Tuning** section to sweep
$C_m$ values and find the one that minimizes NRMSE against COMSOL FSI results.
For the 2-ray fin at $\\Phi = 40°$, the optimal value is approximately $C_m \\approx 0.75$.

**Important:** If $C_n$ was previously calibrated to match CFD with the RFT-only model,
it already absorbs all force physics into one coefficient. For Morison, consider jointly
tuning both $C_n$ and $C_m$.
""")

    st.divider()

    # ─────────── 6. Lighthill ───────────
    st.subheader("6. Model III — Lighthill's Large-Amplitude EBT")
    st.markdown("""
Lighthill's (1971) elongated body theory computes **reactive** forces from the
rate of change of virtual fluid momentum. Unlike RFT/Morison which use local
velocity/acceleration, Lighthill's theory captures the momentum exchange between
the body and the fluid along its entire length.
""")

    st.markdown("**Eq. (L1)** — Reactive force per panel:")
    st.latex(r"""
\mathbf{F}_{\text{reactive}} =
-\frac{d}{dt}\int_{h_{\min}}^{h_{\max}} m(h)\, v_n(h)\, \mathbf{e}_n(h)\; dh
\;+\; \tfrac{1}{2}\, m(h_{\max})\, v_n^2(h_{\max})\, \mathbf{e}_n(h_{\max})
""")

    col_l1, col_l2 = st.columns(2)
    with col_l1:
        st.markdown("**Term 1 — Bulk momentum rate:**")
        st.markdown("""
The time derivative of the total virtual fluid momentum
carried by the membrane. This is the dominant reactive force.
""")
    with col_l2:
        st.markdown("**Term 2 — Trailing-edge shedding:**")
        st.markdown("""
Momentum shed into the wake at the trailing edge
($h = h_{\\max}$). This represents the propulsive jet.
""")

    st.markdown("**Eq. (L2)** — Added mass per unit height — flat plate model:")
    st.latex(r"""
m(h) = \frac{\pi}{4}\,\rho\, W(h)^2
""")
    st.caption("where $W(h) = w_{\\max}(h) = \\sqrt{2h^2(1-\\cos\\Delta q)+D^2}$ is the local chord.")

    st.markdown("**Eq. (L3)** — Momentum rate approximation:")
    st.latex(r"""
\frac{d}{dt}\bigl[m(h)\,v_n\,\mathbf{e}_n\bigr]
\approx m(h)\,\dot{v}_n\,\mathbf{e}_n
""")
    st.caption(
        "Higher-order terms from $\\dot{m}$ and $\\dot{\\mathbf{e}}_n$ "
        "are small for the quasi-steady regime."
    )

    with st.expander("Key differences: RFT vs Morison vs Lighthill"):
        st.markdown("""
| Aspect | RFT | Morison | Lighthill |
|--------|-----|---------|-----------|
| **Force type** | Resistive (drag) | Resistive + reactive | Reactive |
| **Depends on** | $v_n^2$ | $v_n^2$ and $\\dot{v}_n$ | $v_n$, $\\dot{v}_n$, geometry |
| **Physical model** | Flat plate drag | Flat plate drag + added mass | Elongated body momentum |
| **Parameters** | $C_n$ | $C_n$, $C_m$ | $\\rho$ only (no empirical coefficients) |
| **Best for** | Low frequency, small $\\Delta q$ | Moderate frequency | Large amplitude, high frequency |
| **Limitations** | No inertia effects | Empirical coefficients needed | Slender body assumption |
| **Propulsive force** | From drag asymmetry | From drag + inertia asymmetry | From trailing-edge shedding |
""")

    st.divider()

    # ─────────── 7. Force Components ───────────
    st.subheader("7. Force Component Mapping")
    st.markdown("""
The three force components correspond to physical directions in the robot frame:
""")
    st.latex(r"""
\mathbf{F} =
\begin{bmatrix} F_{\text{heave}} \\ F_{\text{sway}} \\ F_{\text{surge}} \end{bmatrix}
=
\begin{bmatrix}
\text{Vertical (lift)} & \leftarrow \sin(q) \text{ terms} \\
\text{Lateral} & \leftarrow \cos(q) \text{ terms} \\
\text{Propulsive (thrust)} & \leftarrow \sin(\Delta q) \text{ terms}
\end{bmatrix}
""")

    st.markdown("""
| Component | Python index | Physical meaning | Dominant term |
|-----------|:-----------:|------------------|---------------|
| **Heave** $F_x$ | `F[0]` | Vertical lift | $\\sin q_i$ |
| **Sway** $F_y$ | `F[1]` | Lateral force | $\\cos q_i$ |
| **Surge** $F_z$ | `F[2]` | Propulsive thrust | $\\sin \\Delta q$ |
""")

    st.divider()

    # ─────────── 8. References ───────────
    st.subheader("8. References")
    st.markdown("""
1. **Ali, A.I.** (2021). *Kinematic and Dynamic Modeling of Biomimetic Undulating
   Fins for Underwater Robots*. MSc Thesis, Bahçeşehir University.

2. **Morison, J.R., Johnson, J.W., Schaaf, S.A.** (1950). *The force exerted by
   surface waves on piles*. Journal of Petroleum Technology, 2(05), 149–154.

3. **Lighthill, M.J.** (1971). *Large-amplitude elongated-body theory of fish
   locomotion*. Proceedings of the Royal Society B, 179(1055), 125–138.

4. **Taylor, G.I.** (1952). *Analysis of the swimming of long and narrow animals*.
   Proceedings of the Royal Society A, 214(1117), 158–183.

5. **Sarpkaya, T.** (2010). *Wave Forces on Offshore Structures*. Cambridge
   University Press.

6. **Lighthill, M.J.** (1960). *Note on the swimming of slender fish*. Journal
   of Fluid Mechanics, 9(2), 305–317.
""")


# ──────────────────────────── main ───────────────────────────────────

def main():
    st.title("🐟 Undulating Fin Dynamics")
    st.caption("Interactive kinematic & dynamic modeling of biomimetic undulating fins")

    params, sim, hydro_method = build_sidebar()

    tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9 = st.tabs([
        "📈 Time Domain",
        "🔄 Parametric Sweep",
        "🤖 Robot Simulation",
        "🌊 Kinematics",
        "📊 CFD Comparison",
        "🔬 OpenFOAM",
        "📐 Math Models",
        "🐛 MATLAB Corrections",
        "ℹ️ About",
    ])

    with tab1:
        tab_time_domain(params, sim, hydro_method)
    with tab2:
        tab_parametric(params, sim, hydro_method)
    with tab3:
        tab_robot(params, sim)
    with tab4:
        tab_kinematics(params, sim)
    with tab5:
        tab_cfd_comparison()
    with tab6:
        tab_openfoam(params)
    with tab7:
        tab_math_models()
    with tab8:
        tab_matlab_corrections()
    with tab9:
        tab_about()


if __name__ == "__main__":
    main()
