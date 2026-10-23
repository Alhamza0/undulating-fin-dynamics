#!/usr/bin/env python3
"""
Generate two scientific diagrams for the undulating fin paper:
  1. Fin schematic – rays, membrane, kinematic variables
  2. Force physics – RFT drag vs Morison inertia, thrust mechanism
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Arc
from pathlib import Path

SRC = Path(__file__).parent / "src"
sys.path.insert(0, str(SRC))

# ── Output directory ────────────────────────────────────────────────
FIG_DIR = Path(__file__).parent / "paper" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Style (match existing paper figures) ────────────────────────────
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# ════════════════════════════════════════════════════════════════════
# FIGURE 1: Fin Schematic Diagram
# ════════════════════════════════════════════════════════════════════

def fig_fin_schematic():
    """
    Side-view schematic of the undulating fin showing:
    - Body axis with pivot joints
    - N rays at different phase angles (traveling wave)
    - Membrane between ray tips
    - Annotations: D, h_max, A, Φ, wave direction
    - Coordinate system (Heave, Sway, Surge)
    """
    fig, ax = plt.subplots(figsize=(7.5, 4.5))

    # ── Parameters ──────────────────────────────────────────────────
    N = 5                       # number of rays to show
    D = 1.0                     # normalized spacing
    L_ray = 1.5                 # ray length (visual)
    A_deg = 30.0                # amplitude in degrees
    Phi_deg = 40.0              # inter-ray phase lag
    t = 0.15                    # time instant (chosen for nice wave shape)
    f = 1.0                     # frequency

    A = np.radians(A_deg)
    Phi = np.radians(Phi_deg)

    # Compute ray angles at this instant
    angles = [A * np.sin(2 * np.pi * f * t - i * Phi) for i in range(N)]

    # Pivot positions along body axis
    pivots_x = [i * D for i in range(N)]
    pivots_y = [0.0] * N

    # Ray tip positions
    tips_x = [pivots_x[i] + L_ray * np.sin(angles[i]) for i in range(N)]
    tips_y = [pivots_y[i] + L_ray * np.cos(angles[i]) for i in range(N)]

    # ── Body axis (base plate / fin body) ───────────────────────────
    body_xmin = -0.6
    body_xmax = pivots_x[-1] + 0.6
    ax.plot([body_xmin, body_xmax], [0, 0], 'k-', lw=3, zorder=5)
    ax.text((body_xmin + body_xmax) / 2, -0.23,
            "Fin body axis", ha="center", fontsize=10, fontstyle="italic",
            color="0.3")

    # ── Draw rays and pivots ────────────────────────────────────────
    ray_color = "#1565C0"
    for i in range(N):
        # Pivot joint (circle)
        ax.plot(pivots_x[i], pivots_y[i], 'o', ms=7,
                mfc="white", mec="black", mew=1.5, zorder=10)
        # Ray (thick line)
        ax.plot([pivots_x[i], tips_x[i]], [pivots_y[i], tips_y[i]],
                '-', lw=2.5, color=ray_color, zorder=8, solid_capstyle='round')
        # Ray tip
        ax.plot(tips_x[i], tips_y[i], 'o', ms=4, mfc=ray_color,
                mec=ray_color, zorder=9)
        # Label ray
        ax.text(pivots_x[i], -0.12, f"$r_{i+1}$", ha="center",
                fontsize=9, color="0.3")

    # ── Draw membrane (filled region between tips) ──────────────────
    membrane_color = "#90CAF9"
    for i in range(N - 1):
        # Create membrane as a filled polygon between two rays
        n_pts = 20
        # Interpolate along each ray from pivot to tip
        for frac in np.linspace(0.3, 1.0, 8):
            x1 = pivots_x[i] + frac * L_ray * np.sin(angles[i])
            y1 = pivots_y[i] + frac * L_ray * np.cos(angles[i])
            x2 = pivots_x[i+1] + frac * L_ray * np.sin(angles[i+1])
            y2 = pivots_y[i+1] + frac * L_ray * np.cos(angles[i+1])
            ax.plot([x1, x2], [y1, y2], '-', lw=0.5, color=membrane_color,
                    alpha=0.6, zorder=3)

        # Fill membrane polygon
        xs = []
        ys = []
        fracs = np.linspace(0.25, 1.0, 15)
        for frac in fracs:
            xs.append(pivots_x[i] + frac * L_ray * np.sin(angles[i]))
            ys.append(pivots_y[i] + frac * L_ray * np.cos(angles[i]))
        for frac in fracs[::-1]:
            xs.append(pivots_x[i+1] + frac * L_ray * np.sin(angles[i+1]))
            ys.append(pivots_y[i+1] + frac * L_ray * np.cos(angles[i+1]))
        ax.fill(xs, ys, color=membrane_color, alpha=0.25, zorder=2)

    # Draw tip-to-tip curve (envelope)
    from scipy.interpolate import CubicSpline
    cs = CubicSpline(pivots_x, tips_x)
    cs_y = CubicSpline(pivots_x, tips_y)
    x_fine = np.linspace(pivots_x[0], pivots_x[-1], 200)
    ax.plot(cs(x_fine), cs_y(x_fine), '--', color=ray_color, lw=1.2,
            alpha=0.7, zorder=7, label="Ray-tip envelope")

    # ── Annotation: D (ray spacing) ────────────────────────────────
    y_D = -0.38
    i_D = 1  # between ray 2 and ray 3
    ax.annotate("", xy=(pivots_x[i_D+1], y_D), xytext=(pivots_x[i_D], y_D),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text((pivots_x[i_D] + pivots_x[i_D+1]) / 2, y_D - 0.09,
            "$D$", ha="center", fontsize=12, fontweight="bold")

    # ── Annotation: h_max (membrane height) ─────────────────────────
    i_h = 2
    # Draw along the ray
    hx0 = pivots_x[i_h] + 0.25 * L_ray * np.sin(angles[i_h])
    hy0 = pivots_y[i_h] + 0.25 * L_ray * np.cos(angles[i_h])
    hx1 = tips_x[i_h]
    hy1 = tips_y[i_h]
    # Arrow along ray for h_max
    offset = 0.15
    dx_perp = -np.cos(angles[i_h]) * offset
    dy_perp = np.sin(angles[i_h]) * offset
    ax.annotate("", xy=(hx1 + dx_perp, hy1 + dy_perp),
                xytext=(hx0 + dx_perp, hy0 + dy_perp),
                arrowprops=dict(arrowstyle="<->", color="#D32F2F", lw=1.2))
    mid_hx = (hx0 + hx1) / 2 + dx_perp - 0.12
    mid_hy = (hy0 + hy1) / 2 + dy_perp
    ax.text(mid_hx, mid_hy, "$h_{\\mathrm{max}}$", fontsize=11,
            color="#D32F2F", fontweight="bold", ha="center")

    # ── Annotation: A (amplitude arc) ───────────────────────────────
    i_A = 0
    arc_r = 0.7
    # Draw arc from 0 to angles[i_A]
    theta_start = 90 - np.degrees(A)
    theta_end = 90 + np.degrees(A)
    arc = Arc((pivots_x[i_A], pivots_y[i_A]), 2*arc_r, 2*arc_r,
              angle=0, theta1=theta_start, theta2=theta_end,
              color="#388E3C", lw=1.5, linestyle="--")
    ax.add_patch(arc)
    # Ghost ray at max angle
    ax.plot([pivots_x[i_A], pivots_x[i_A] + L_ray * np.sin(A)],
            [pivots_y[i_A], pivots_y[i_A] + L_ray * np.cos(A)],
            ':', lw=1.2, color="#388E3C", alpha=0.5, zorder=4)
    ax.plot([pivots_x[i_A], pivots_x[i_A] + L_ray * np.sin(-A)],
            [pivots_y[i_A], pivots_y[i_A] + L_ray * np.cos(-A)],
            ':', lw=1.2, color="#388E3C", alpha=0.5, zorder=4)
    # Vertical reference line
    ax.plot([pivots_x[i_A], pivots_x[i_A]],
            [pivots_y[i_A], pivots_y[i_A] + arc_r + 0.1],
            ':', lw=0.8, color="gray", alpha=0.5, zorder=4)
    # Label A
    ax.text(pivots_x[i_A] - 0.35, 0.65, "$A$", fontsize=12,
            color="#388E3C", fontweight="bold")

    # ── Annotation: Φ (phase between ray 1 and ray 2) ──────────────
    i_P = 0
    phi_r = 0.5
    # Show angle difference between ray 1 and ray 2
    th1 = 90 - np.degrees(angles[i_P])
    th2 = 90 - np.degrees(angles[i_P + 1])
    theta_lo = min(th1, th2)
    theta_hi = max(th1, th2)
    arc_phi = Arc((pivots_x[i_P + 1], pivots_y[i_P + 1]),
                  2 * phi_r, 2 * phi_r,
                  angle=0, theta1=theta_lo, theta2=theta_hi,
                  color="#E65100", lw=1.5)
    ax.add_patch(arc_phi)
    ax.text(pivots_x[i_P + 1] + 0.08, 0.55, "$\\Phi$", fontsize=12,
            color="#E65100", fontweight="bold")

    # ── Wave propagation arrow ──────────────────────────────────────
    wave_y = 1.85
    ax.annotate("", xy=(pivots_x[-1] + 0.3, wave_y),
                xytext=(pivots_x[0] - 0.1, wave_y),
                arrowprops=dict(arrowstyle="->,head_width=0.15,head_length=0.15",
                                color="#6A1B9A", lw=2))
    ax.text((pivots_x[0] + pivots_x[-1]) / 2, wave_y + 0.12,
            "Wave propagation", ha="center", fontsize=10,
            color="#6A1B9A", fontstyle="italic")

    # ── Coordinate system (inset) ───────────────────────────────────
    cx, cy = body_xmax + 0.1, 1.0
    arr_len = 0.45
    head = dict(arrowstyle="->,head_width=0.12,head_length=0.1", lw=1.5)
    # Heave (up)
    ax.annotate("", xy=(cx, cy + arr_len), xytext=(cx, cy),
                arrowprops=dict(**head, color="#D32F2F"))
    ax.text(cx + 0.08, cy + arr_len + 0.02, "Heave", fontsize=8,
            color="#D32F2F", va="bottom")
    # Sway (right along body)
    ax.annotate("", xy=(cx + arr_len, cy), xytext=(cx, cy),
                arrowprops=dict(**head, color="#1565C0"))
    ax.text(cx + arr_len + 0.05, cy - 0.02, "Sway", fontsize=8,
            color="#1565C0", va="top")
    # Surge (out of page, circle with dot)
    ax.plot(cx, cy, 'o', ms=10, mfc="white", mec="#388E3C", mew=1.5, zorder=15)
    ax.plot(cx, cy, '.', ms=4, color="#388E3C", zorder=16)
    ax.text(cx - 0.3, cy + 0.12, "Surge\n(thrust)", fontsize=8,
            color="#388E3C", ha="center")

    # ── Final formatting ────────────────────────────────────────────
    ax.set_xlim(body_xmin - 0.5, body_xmax + 1.1)
    ax.set_ylim(-0.6, 2.15)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fin_schematic.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fin_schematic.png", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ fin_schematic.pdf / .png")


# ════════════════════════════════════════════════════════════════════
# FIGURE 2: Force Physics Illustration
# ════════════════════════════════════════════════════════════════════

def fig_force_physics():
    """
    Two-panel illustration comparing RFT drag and Morison inertia:
    (a) Membrane element with v_n, RFT drag f_n opposing motion
    (b) Membrane element with dv_n/dt, Morison inertia force, thrust arrow
    Plus a bottom panel showing force decomposition into Heave/Sway/Surge.
    """
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 4.2))

    for ax in axes:
        ax.set_aspect("equal")
        ax.axis("off")

    # ── Common geometry: two adjacent rays with membrane ────────────
    def draw_panel(ax, panel_label, title):
        """Draw two rays with membrane panel."""
        # Ray geometry
        px1, py1 = 0.0, 0.0      # pivot 1
        px2, py2 = 1.3, 0.0      # pivot 2
        L = 2.0                    # ray length
        q1 = np.radians(25)      # ray 1 angle
        q2 = np.radians(-10)     # ray 2 angle

        # Tip positions
        tx1 = px1 + L * np.sin(q1)
        ty1 = py1 + L * np.cos(q1)
        tx2 = px2 + L * np.sin(q2)
        ty2 = py2 + L * np.cos(q2)

        # Body axis
        ax.plot([-0.4, 2.0], [0, 0], 'k-', lw=3, zorder=5)

        # Pivot circles
        for px, py in [(px1, py1), (px2, py2)]:
            ax.plot(px, py, 'o', ms=6, mfc="white", mec="black",
                    mew=1.3, zorder=10)

        # Rays
        ray_c = "#1565C0"
        ax.plot([px1, tx1], [py1, ty1], '-', lw=2.5, color=ray_c, zorder=8)
        ax.plot([px2, tx2], [py2, ty2], '-', lw=2.5, color=ray_c, zorder=8)
        ax.plot(tx1, ty1, 'o', ms=3, color=ray_c, zorder=9)
        ax.plot(tx2, ty2, 'o', ms=3, color=ray_c, zorder=9)

        # Membrane fill
        fracs = np.linspace(0.25, 1.0, 20)
        xs_left = [px1 + f * L * np.sin(q1) for f in fracs]
        ys_left = [py1 + f * L * np.cos(q1) for f in fracs]
        xs_right = [px2 + f * L * np.sin(q2) for f in fracs[::-1]]
        ys_right = [py2 + f * L * np.cos(q2) for f in fracs[::-1]]
        ax.fill(xs_left + xs_right, ys_left + ys_right,
                color="#90CAF9", alpha=0.25, zorder=2)

        # Membrane grid lines
        for frac in np.linspace(0.3, 1.0, 6):
            x1 = px1 + frac * L * np.sin(q1)
            y1 = py1 + frac * L * np.cos(q1)
            x2 = px2 + frac * L * np.sin(q2)
            y2 = py2 + frac * L * np.cos(q2)
            ax.plot([x1, x2], [y1, y2], '-', lw=0.5, color="#90CAF9",
                    alpha=0.5, zorder=3)

        # Labels
        ax.text(px1, -0.17, "$r_i$", ha="center", fontsize=9, color="0.3")
        ax.text(px2, -0.17, "$r_{i+1}$", ha="center", fontsize=9, color="0.3")

        # Title
        ax.set_title(title, fontsize=11, fontweight="bold", pad=10)

        # Panel label
        ax.text(-0.45, 2.35, panel_label, fontsize=13, fontweight="bold",
                va="top")

        return q1, q2, px1, py1, px2, py2, tx1, ty1, tx2, ty2, L

    # ════════════════════════════════════════════════════════════════
    # Panel (a): RFT – Drag Force
    # ════════════════════════════════════════════════════════════════
    ax = axes[0]
    q1, q2, px1, py1, px2, py2, tx1, ty1, tx2, ty2, L = \
        draw_panel(ax, "(a)", "RFT — Resistive Drag")

    # Midpoint on membrane
    frac_pt = 0.65
    mx = (px1 + frac_pt * L * np.sin(q1) + px2 + frac_pt * L * np.sin(q2)) / 2
    my = (py1 + frac_pt * L * np.cos(q1) + py2 + frac_pt * L * np.cos(q2)) / 2

    # Normal direction to membrane at this point
    # Membrane tangent ≈ direction from ray_i to ray_{i+1} at this height
    tang_x = (px2 + frac_pt * L * np.sin(q2)) - (px1 + frac_pt * L * np.sin(q1))
    tang_y = (py2 + frac_pt * L * np.cos(q2)) - (py1 + frac_pt * L * np.cos(q1))
    tang_len = np.sqrt(tang_x**2 + tang_y**2)
    tang_x /= tang_len
    tang_y /= tang_len
    # Normal = perpendicular to tangent
    norm_x = -tang_y
    norm_y = tang_x

    # Velocity arrow (membrane moves in normal direction)
    v_len = 0.6
    arr_kw = dict(arrowstyle="->,head_width=0.1,head_length=0.08", lw=1.8)
    ax.annotate("", xy=(mx + norm_x * v_len, my + norm_y * v_len),
                xytext=(mx, my),
                arrowprops=dict(**arr_kw, color="#7B1FA2"))
    ax.text(mx + norm_x * v_len + 0.08, my + norm_y * v_len + 0.08,
            "$v_n$", fontsize=11, color="#7B1FA2", fontweight="bold")

    # RFT drag force (opposes velocity)
    f_len = 0.5
    ax.annotate("", xy=(mx - norm_x * f_len, my - norm_y * f_len),
                xytext=(mx, my),
                arrowprops=dict(**arr_kw, color="#D32F2F"))
    ax.text(mx - norm_x * f_len - 0.05, my - norm_y * f_len - 0.15,
            "$f_n = -\\frac{1}{2}\\rho C_n v_n^2$",
            fontsize=9, color="#D32F2F", fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="#D32F2F",
                      boxstyle="round,pad=0.2", alpha=0.9))

    # Small element highlight
    elem_w = 0.12
    rect_pts_x = [mx - tang_x * elem_w - norm_x * 0.03,
                  mx + tang_x * elem_w - norm_x * 0.03,
                  mx + tang_x * elem_w + norm_x * 0.03,
                  mx - tang_x * elem_w + norm_x * 0.03]
    rect_pts_y = [my - tang_y * elem_w - norm_y * 0.03,
                  my + tang_y * elem_w - norm_y * 0.03,
                  my + tang_y * elem_w + norm_y * 0.03,
                  my - tang_y * elem_w + norm_y * 0.03]
    ax.fill(rect_pts_x, rect_pts_y, color="#FFCDD2", alpha=0.7, zorder=6)
    ax.plot(rect_pts_x + [rect_pts_x[0]], rect_pts_y + [rect_pts_y[0]],
            '-', color="#D32F2F", lw=1, zorder=6)
    ax.text(mx + tang_x * 0.34, my + tang_y * 0.34 + 0.06,
            "$dA$", fontsize=8, color="#D32F2F", fontstyle="italic",
            ha="left", va="bottom", zorder=7)

    # Drag alone already produces net thrust once the rays are phase
    # shifted -- the travelling wave is what breaks time-reversal
    # symmetry, not the inertia term.  Panel (a) therefore carries the
    # same thrust arrow as panel (b); omitting it here was what made
    # this figure wrongly imply that thrust needs the inertia term.
    thrust_x_a = mx + 0.62
    thrust_y_a = my - 1.05
    ax.annotate("",
                xy=(thrust_x_a + 0.5, thrust_y_a),
                xytext=(thrust_x_a, thrust_y_a),
                arrowprops=dict(arrowstyle="->,head_width=0.15,head_length=0.12",
                                color="#388E3C", lw=2.5))
    ax.text(thrust_x_a + 0.55, thrust_y_a + 0.02,
            "Net\nthrust", fontsize=9, color="#388E3C",
            fontweight="bold", va="center")

    # Key point text
    ax.text(0.65, -1.00,
            "Drag: in phase with $v_n$\n"
            "Sets the mean thrust ($\\propto C_n$)",
            fontsize=8, ha="center", color="0.3",
            fontstyle="italic",
            bbox=dict(facecolor="#FFF9C4", edgecolor="0.7",
                      boxstyle="round,pad=0.3", alpha=0.8))

    ax.set_xlim(-0.55, 2.35)
    ax.set_ylim(-1.30, 2.5)

    # ════════════════════════════════════════════════════════════════
    # Panel (b): Morison – Drag + Inertia
    # ════════════════════════════════════════════════════════════════
    ax = axes[1]
    q1, q2, px1, py1, px2, py2, tx1, ty1, tx2, ty2, L = \
        draw_panel(ax, "(b)", "Morison — Drag + Inertia")

    # Same midpoint
    mx = (px1 + frac_pt * L * np.sin(q1) + px2 + frac_pt * L * np.sin(q2)) / 2
    my = (py1 + frac_pt * L * np.cos(q1) + py2 + frac_pt * L * np.cos(q2)) / 2

    # Recompute normal
    tang_x = (px2 + frac_pt * L * np.sin(q2)) - (px1 + frac_pt * L * np.sin(q1))
    tang_y = (py2 + frac_pt * L * np.cos(q2)) - (py1 + frac_pt * L * np.cos(q1))
    tang_len = np.sqrt(tang_x**2 + tang_y**2)
    tang_x /= tang_len
    tang_y /= tang_len
    norm_x = -tang_y
    norm_y = tang_x

    # Small element highlight (recompute for this panel)
    rect_pts_x2 = [mx - tang_x * elem_w - norm_x * 0.03,
                   mx + tang_x * elem_w - norm_x * 0.03,
                   mx + tang_x * elem_w + norm_x * 0.03,
                   mx - tang_x * elem_w + norm_x * 0.03]
    rect_pts_y2 = [my - tang_y * elem_w - norm_y * 0.03,
                   my + tang_y * elem_w - norm_y * 0.03,
                   my + tang_y * elem_w + norm_y * 0.03,
                   my - tang_y * elem_w + norm_y * 0.03]
    ax.fill(rect_pts_x2, rect_pts_y2, color="#FFCDD2", alpha=0.7, zorder=6)
    ax.plot(rect_pts_x2 + [rect_pts_x2[0]], rect_pts_y2 + [rect_pts_y2[0]],
            '-', color="#D32F2F", lw=1, zorder=6)

    # Velocity arrow
    ax.annotate("", xy=(mx + norm_x * v_len, my + norm_y * v_len),
                xytext=(mx, my),
                arrowprops=dict(**arr_kw, color="#7B1FA2"))
    ax.text(mx + norm_x * v_len + 0.08, my + norm_y * v_len + 0.08,
            "$v_n$", fontsize=11, color="#7B1FA2", fontweight="bold")

    # Drag component (same as RFT, slightly shorter)
    f_drag_len = 0.35
    drag_end_x = mx - norm_x * f_drag_len
    drag_end_y = my - norm_y * f_drag_len
    ax.annotate("", xy=(drag_end_x, drag_end_y), xytext=(mx, my),
                arrowprops=dict(**arr_kw, color="#D32F2F"))

    # Inertia component (shifted along normal, different color)
    # Acceleration arrow (dv_n/dt)
    acc_x = norm_x * 0.4
    acc_y = norm_y * 0.4
    acc_start_x = mx + tang_x * 0.25
    acc_start_y = my + tang_y * 0.25
    ax.annotate("",
                xy=(acc_start_x + acc_x, acc_start_y + acc_y),
                xytext=(acc_start_x, acc_start_y),
                arrowprops=dict(**arr_kw, color="#FF6F00"))
    ax.text(acc_start_x + acc_x + 0.05, acc_start_y + acc_y + 0.05,
            "$\\dot{v}_n$", fontsize=10, color="#FF6F00", fontweight="bold")

    # Inertia force (opposes acceleration)
    inertia_len = 0.35
    inertia_end_x = acc_start_x - norm_x * inertia_len
    inertia_end_y = acc_start_y - norm_y * inertia_len
    ax.annotate("",
                xy=(inertia_end_x, inertia_end_y),
                xytext=(acc_start_x, acc_start_y),
                arrowprops=dict(**arr_kw, color="#FF6F00"))

    # Formula box
    ax.text(0.65, -0.62,
            "$F = F_{\\mathrm{drag}} + F_{\\mathrm{inertia}}$\n"
            "drag: $-\\frac{1}{2}\\rho C_n v_n^2$    "
            "inertia: $-\\rho C_m \\frac{\\pi}{4} w \\dot{v}_n$",
            fontsize=8, ha="center",
            bbox=dict(facecolor="white", edgecolor="#FF6F00",
                      boxstyle="round,pad=0.3", alpha=0.9))

    # Thrust arrow (surge component)
    # Show decomposition: the surge component of the inertia force generates net thrust
    thrust_x = mx + 0.62
    thrust_y = my - 1.05
    ax.annotate("",
                xy=(thrust_x + 0.5, thrust_y),
                xytext=(thrust_x, thrust_y),
                arrowprops=dict(arrowstyle="->,head_width=0.15,head_length=0.12",
                                color="#388E3C", lw=2.5))
    ax.text(thrust_x + 0.55, thrust_y + 0.02,
            "Net\nthrust", fontsize=9, color="#388E3C",
            fontweight="bold", va="center")

    # Legend labels
    ax.plot([], [], '-', color="#D32F2F", lw=2, label="Drag ($f_n$)")
    ax.plot([], [], '-', color="#FF6F00", lw=2, label="Inertia ($f_{\\mathrm{inertia}}$)")
    ax.plot([], [], '-', color="#388E3C", lw=2.5, label="Surge (thrust)")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # Key point text
    ax.text(0.65, -1.00,
            "Inertia: in quadrature with drag\n"
            "Sets the waveform phase, not the mean",
            fontsize=8, ha="center", color="0.3", fontstyle="italic",
            bbox=dict(facecolor="#C8E6C9", edgecolor="0.7",
                      boxstyle="round,pad=0.3", alpha=0.8))

    ax.set_xlim(-0.55, 2.35)
    ax.set_ylim(-1.30, 2.5)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "force_physics.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "force_physics.png", bbox_inches="tight")
    plt.close(fig)
    print("  ✓ force_physics.pdf / .png")


# ════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    print("Generating scientific diagrams...")
    fig_fin_schematic()
    fig_force_physics()
    print("Done — figures saved to", FIG_DIR)
