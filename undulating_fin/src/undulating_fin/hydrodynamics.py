"""
Hydrodynamic force and torque computation.

Supports three models:
- **RFT** (Resistive Force Theory): drag-only, quasi-steady (original).
- **Morison**: RFT drag + added-mass inertia term.
- **Lighthill**: Large-amplitude elongated body theory (reactive forces).

References
----------
- Thesis Eq. (12): Drag force per unit area f_n
- Thesis Eq. (13): Total force F_N by double integration over membrane
- Thesis Eq. (14): Normal velocity v_n(h, w̄)
- Thesis Eq. (15): Geometric factor δ(h, w̄)
- Thesis Eq. (16): Normal vector e_n(h, w̄)
- Thesis Sec. 2.3.2: Hydrodynamic forces derivation
- Morison et al. (1950): F = 0.5*ρ*Cd*A*|vn|*vn + ρ*Cm*V*dvn/dt
- Lighthill (1971): Large-amplitude elongated body theory
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import dblquad
from undulating_fin.params import FinParams


def _compute_vn_en_delta(
    h: float,
    w_bar: float,
    q1: float,
    q2: float,
    dq1: float,
    dq2: float,
    D: float,
) -> tuple[float, np.ndarray, float]:
    """Compute normal velocity, normal vector, and δ for a membrane point.

    Parameters
    ----------
    h : float
        Height along membrane [m].
    w_bar : float
        Normalized lengthwise coordinate ∈ [0, 1]. Thesis: w̄ = w/w_max.
    q1, q2 : float
        Joint angles of adjacent rays [rad].
    dq1, dq2 : float
        Joint angular velocities [rad/s].
    D : float
        Inter-ray distance [m].

    Returns
    -------
    vn : float
        Normal velocity component [m/s]. Eq. (14).
    en : ndarray, shape (3,)
        Unit normal vector components [e_x, e_y, e_z]. Eq. (16).
    delta : float
        Geometric scale factor [m]. Eq. (15).
    """
    dq = q1 - q2  # Δq_{1|2}
    cos_dq = np.cos(dq)
    sin_dq = np.sin(dq)

    s1, s2 = np.sin(q1), np.sin(q2)
    c1, c2 = np.cos(q1), np.cos(q2)

    # Eq. (15): δ = sqrt(2D²·w̄(w̄−1)(1−cos Δq) + h²·sin²Δq + D²)
    delta = np.sqrt(
        2 * D**2 * w_bar * (w_bar - 1) * (1 - cos_dq)
        + h**2 * sin_dq**2
        + D**2
    )

    if abs(delta) < 1e-30:
        return 0.0, np.zeros(3), 1e-30

    # Eq. (14): v_n = (D·h/δ)·(dq1 + (1−cos Δq)(dq1+dq2)w̄²
    #                          + (cos Δq·(dq1+dq2) − 2dq1)·w̄)
    vn = (D * h / delta) * (
        dq1
        + (1 - cos_dq) * (dq1 + dq2) * w_bar**2
        + (cos_dq * (dq1 + dq2) - 2 * dq1) * w_bar
    )

    # Eq. (16): e_n = (D/δ)·[s1(w̄−1) − s2·w̄, −c1(w̄−1) + c2·w̄, −sin(Δq)·h/D]
    en = np.array([
        s1 * (w_bar - 1) - s2 * w_bar,
        -c1 * (w_bar - 1) + c2 * w_bar,
        -sin_dq * h / D,
    ]) * (D / delta)

    return vn, en, delta


# ── Valid method names ──
HYDRO_METHODS = ("rft", "morison", "lighthill")


def _compute_vn_dot(
    h: float,
    w_bar: float,
    q1: float,
    q2: float,
    dq1: float,
    dq2: float,
    ddq1: float,
    ddq2: float,
    D: float,
) -> float:
    r"""Time derivative of normal velocity (dv_n/dt) for Morison inertia term.

    Derived by differentiating Eq. (14) w.r.t. time. The geometric factor
    δ also depends on q (hence time), but to leading order we compute
    the dominant acceleration contribution while holding δ constant
    (consistent with Morison's quasi-steady framework).

    The dominant part is:
        dv_n/dt ≈ (D·h/δ) · [ddq1 + (1−cos Δq)(ddq1+ddq2)w̄²
                   + (cos Δq·(ddq1+ddq2) − 2·ddq1)·w̄
                   + sin Δq·(dq1−dq2)·((dq1+dq2)w̄² + (dq1+dq2−2dq1/(cosΔq))·w̄ ...)]

    For clarity and robustness we use the full analytic derivative below.
    """
    dq = q1 - q2
    ddq = dq1 - dq2
    cos_dq = np.cos(dq)
    sin_dq = np.sin(dq)

    delta = np.sqrt(
        2 * D**2 * w_bar * (w_bar - 1) * (1 - cos_dq)
        + h**2 * sin_dq**2
        + D**2
    )
    if abs(delta) < 1e-30:
        return 0.0

    # --- vn expression pieces (Eq. 14) ---
    A_coeff = dq1
    B_coeff = (1 - cos_dq) * (dq1 + dq2)       # multiplied by w̄²
    C_coeff = cos_dq * (dq1 + dq2) - 2 * dq1    # multiplied by w̄

    # --- Time derivatives of A, B, C coefficients ---
    dA = ddq1
    dB = sin_dq * ddq * (dq1 + dq2) + (1 - cos_dq) * (ddq1 + ddq2)
    dC = -sin_dq * ddq * (dq1 + dq2) + cos_dq * (ddq1 + ddq2) - 2 * ddq1

    # Leading-order: hold D·h/δ constant, differentiate the polynomial
    prefix = D * h / delta
    vn_dot = prefix * (
        dA
        + dB * w_bar**2
        + dC * w_bar
    )

    # --- Correction from d(1/δ)/dt (chain rule) ---
    # dδ²/dt = 2D²·w̄(w̄−1)·sin(Δq)·(dq1−dq2) + 2h²·sin(Δq)·cos(Δq)·(dq1−dq2)
    d_delta_sq = (
        2 * D**2 * w_bar * (w_bar - 1) * sin_dq * ddq
        + 2 * h**2 * sin_dq * cos_dq * ddq
    )
    d_delta = d_delta_sq / (2 * delta) if abs(delta) > 1e-30 else 0.0

    inner = A_coeff + B_coeff * w_bar**2 + C_coeff * w_bar
    vn_dot -= (D * h * d_delta / delta**2) * inner

    return vn_dot


def compute_panel_forces(
    D: float,
    hmin: float,
    hmax: float,
    q1: float,
    q2: float,
    dq1: float,
    dq2: float,
    Cn: float,
    rho: float,
    *,
    method: str = "rft",
    ddq1: float = 0.0,
    ddq2: float = 0.0,
    Cm: float = 1.0,
) -> np.ndarray:
    """Compute hydrodynamic forces on a single membrane panel Ω_{i|i+1}.

    This is the Python equivalent of the MATLAB ``hydrodynamicforces_Func.m``.

    The force on each surface element ds is (Thesis Eq. 12):
        f_n = −(1/2)·ρ·C_d·e_n·v_n²·sgn(v_n)

    Total force is obtained by integrating over the membrane surface
    (Thesis Eq. 13):
        F = ∫∫ f_n · ds

    Integration is performed over:
        h  : hmin → hmax  (height along membrane)
        w  : 0 → w_max(h)  (length along membrane strip)

    where w_max(h) = sqrt(2h²(1 − cos Δq) + D²).

    Note: The integration variable w relates to w̄ by w̄ = w/w_max,
    so the Jacobian |dw̄/dw| = 1/w_max is included in the integrand.

    Parameters
    ----------
    D : float
        Inter-ray distance [m].
    hmin, hmax : float
        Membrane mounting bounds [m].
    q1, q2 : float
        Joint angles of the two rays [rad].
    dq1, dq2 : float
        Joint angular velocities [rad/s].
    Cn : float
        Normal drag coefficient [-].
    rho : float
        Fluid density [kg/m³].
    method : str
        Hydrodynamic model: 'rft', 'morison', or 'lighthill'.
    ddq1, ddq2 : float
        Joint angular accelerations [rad/s²] (needed for morison/lighthill).
    Cm : float
        Added mass coefficient (Morison only).

    Returns
    -------
    F : ndarray, shape (3,)
        Force vector [F_heave, F_sway, F_surge] [N].
        Component 0 (Fx): heave — sin(q) terms
        Component 1 (Fy): sway  — cos(q) terms
        Component 2 (Fz): surge — sin(Δq) terms (propulsive)
    """
    if method == "lighthill":
        return _lighthill_panel_forces(
            D, hmin, hmax, q1, q2, dq1, dq2, ddq1, ddq2, rho,
        )

    dq = q1 - q2
    cos_dq = np.cos(dq)

    def w_max(h: float) -> float:
        """Upper integration limit. Width of membrane strip at height h."""
        return np.sqrt(2 * h**2 * (1 - cos_dq) + D**2)

    # Pre-compute sin/cos for inner function
    sin_dq = np.sin(dq)
    s1, c1 = np.sin(q1), np.cos(q1)
    s2, c2 = np.sin(q2), np.cos(q2)

    use_morison = (method == "morison")

    def _integrand(w: float, h: float, component: int) -> float:
        """Integrand for force component (0=heave, 1=sway, 2=surge).

        Matches MATLAB convention: F(1)=F_heave, F(2)=F_sway, F(3)=F_surge.
        """
        wm = w_max(h)
        if wm < 1e-30:
            return 0.0
        w_bar = w / wm

        # Eq. (15)
        delta = np.sqrt(
            2 * D**2 * w_bar * (w_bar - 1) * (1 - cos_dq)
            + h**2 * sin_dq**2
            + D**2
        )
        if abs(delta) < 1e-30:
            return 0.0

        # Eq. (14)
        vn = (D * h / delta) * (
            dq1
            + (1 - cos_dq) * (dq1 + dq2) * w_bar**2
            + (cos_dq * (dq1 + dq2) - 2 * dq1) * w_bar
        )

        # Eq. (16) — component of normal vector
        en_components = [
            s1 * (w_bar - 1) - s2 * w_bar,           # e_n,x (heave)
            -c1 * (w_bar - 1) + c2 * w_bar,           # e_n,y (sway)
            -sin_dq * h / D,                           # e_n,z (surge)
        ]
        en_k = en_components[component] * (D / delta)

        # Eq. (12): f_n = −(1/2)·ρ·Cn·en·vn²·sgn(vn)  [drag term]
        sgn_vn = np.sign(vn) if abs(vn) > 1e-30 else 0.0
        f_component = -0.5 * rho * Cn * en_k * vn**2 * sgn_vn

        # Morison inertia term (added-mass reaction force per unit area).
        # For a thin plate of chord c oscillating normally in still fluid,
        # the added mass per unit span is  m_a = ρ·(π/4)·c².
        # Per unit area (÷ c) this gives  f_I = -ρ·Cm·(π/4)·c · dvn/dt.
        # The negative sign arises because the body moves (vn) in still
        # fluid, so the relative flow velocity is U = -vn, and the
        # inertia reaction opposes the body's acceleration.
        if use_morison:
            vn_dot = _compute_vn_dot(
                h, w_bar, q1, q2, dq1, dq2, ddq1, ddq2, D,
            )
            f_component -= rho * Cm * (np.pi / 4) * wm * en_k * vn_dot

        return f_component

    F = np.zeros(3)
    for k in range(3):
        val, _ = dblquad(
            lambda w, h: _integrand(w, h, k),
            hmin,
            hmax,
            0,
            w_max,
            epsabs=1e-8,
            epsrel=1e-8,
        )
        F[k] = val

    return F


def compute_total_forces(
    params: FinParams,
    q: np.ndarray,
    dq: np.ndarray,
    *,
    ddq: np.ndarray | None = None,
    method: str = "rft",
) -> np.ndarray:
    """Sum hydrodynamic forces over all (N-1) membrane panels.

    Thesis Eq. (13):
        F_N = Sigma_{i=1}^{N-1} double_integral f_n ds  over Omega_{i|i+1}

    This replaces the hardcoded 11-call pattern in the MATLAB code.

    Parameters
    ----------
    params : FinParams
        Fin parameters.
    q : 1-D array, shape (N,)
        Joint angles at current time step [rad].
    dq : 1-D array, shape (N,)
        Joint velocities at current time step [rad/s].
    ddq : 1-D array, shape (N,), optional
        Joint accelerations [rad/s^2]. Required for morison/lighthill.
    method : str
        Hydrodynamic model: 'rft', 'morison', or 'lighthill'.

    Returns
    -------
    F_total : ndarray, shape (3,)
        Total force [F_heave, F_sway, F_surge] [N].
    """
    N = len(q)
    if ddq is None:
        if method != "rft":
            # Silently substituting zeros drops the added-mass term, so
            # 'morison' degenerates into 'rft' and reports a plausible
            # but wrong answer — the paper's central result is that the
            # two differ.  Refuse instead.
            raise ValueError(
                f"method={method!r} needs joint accelerations; pass "
                "ddq=... (compute_kinematics returns it as its third "
                "value). Only 'rft' can be evaluated without them."
            )
        ddq = np.zeros(N)
    F_total = np.zeros(3)
    for i in range(N - 1):
        F_panel = compute_panel_forces(
            params.D, params.hmin, params.hmax,
            q[i], q[i + 1], dq[i], dq[i + 1],
            params.Cn, params.rho,
            method=method,
            ddq1=ddq[i], ddq2=ddq[i + 1],
            Cm=params.Cm,
        )
        F_total += F_panel
    return F_total


# ──────────────────── Lighthill's Large-Amplitude EBT ────────────────
def _lighthill_panel_forces(
    D: float,
    hmin: float,
    hmax: float,
    q1: float,
    q2: float,
    dq1: float,
    dq2: float,
    ddq1: float,
    ddq2: float,
    rho: float,
) -> np.ndarray:
    r"""Lighthill's large-amplitude elongated body theory for one panel.

    The reactive force arises from the rate of change of the virtual
    (added-mass) fluid momentum. For an undulating fin membrane the
    added mass per unit height is modelled as that of a flat plate:

        m(h) = (pi/4) * rho * W(h)^2

    where W(h) is the local chord (w_max). The reactive force is:

        F_reactive = - d/dt integral_{hmin}^{hmax} m(h) * v_n(h) * e_n(h) dh
                     + 0.5 * m(hmax) * v_n(hmax)^2 * e_n(hmax)

    The first term is the bulk momentum rate, and the second is the
    trailing-edge momentum shedding term (Lighthill 1971).

    The full time derivative uses the product rule:
      d/dt [m*vn*en] = (dm/dt)*vn*en + m*(dvn/dt)*en + m*vn*(den/dt)
    All three terms are computed analytically.
    """
    dq = q1 - q2
    ddq_diff = dq1 - dq2          # d(Δq)/dt
    cos_dq = np.cos(dq)
    sin_dq = np.sin(dq)
    s1, c1 = np.sin(q1), np.cos(q1)
    s2, c2 = np.sin(q2), np.cos(q2)

    def w_max(h):
        return np.sqrt(2 * h**2 * (1 - cos_dq) + D**2)

    def _added_mass_per_h(h):
        """Added mass per unit height for a flat plate of chord W(h)."""
        W = w_max(h)
        return 0.25 * np.pi * rho * W**2

    def _dm_dt_per_h(h):
        """Time derivative of added mass: dm/dt = (pi/2)*rho*W*dW/dt."""
        W = w_max(h)
        if abs(W) < 1e-30:
            return 0.0
        # dW/dt = h^2 * sin(Δq) * d(Δq)/dt / W
        dW_dt = h**2 * sin_dq * ddq_diff / W
        return 0.5 * np.pi * rho * W * dW_dt

    def _vn_en_at(h, w_bar_eval):
        """Compute vn and en at a specific (h, w_bar_eval) point."""
        delta = np.sqrt(
            2 * D**2 * w_bar_eval * (w_bar_eval - 1) * (1 - cos_dq)
            + h**2 * sin_dq**2 + D**2
        )
        if abs(delta) < 1e-30:
            return 0.0, np.zeros(3)
        vn = (D * h / delta) * (
            dq1
            + (1 - cos_dq) * (dq1 + dq2) * w_bar_eval**2
            + (cos_dq * (dq1 + dq2) - 2 * dq1) * w_bar_eval
        )
        en = np.array([
            s1 * (w_bar_eval - 1) - s2 * w_bar_eval,
            -c1 * (w_bar_eval - 1) + c2 * w_bar_eval,
            -sin_dq * h / D,
        ]) * (D / delta)
        return vn, en

    def _den_dt_at(h, w_bar_eval):
        """Time derivative of normal vector den/dt (analytic)."""
        delta = np.sqrt(
            2 * D**2 * w_bar_eval * (w_bar_eval - 1) * (1 - cos_dq)
            + h**2 * sin_dq**2 + D**2
        )
        if abs(delta) < 1e-30:
            return np.zeros(3)

        # en_unnorm = D * [s1(w-1)-s2*w, -c1(w-1)+c2*w, -sin(dq)*h/D]
        # den_unnorm/dt contributions from dq1, dq2 (angular velocities)
        den_unnorm = D * np.array([
            c1 * dq1 * (w_bar_eval - 1) - c2 * dq2 * w_bar_eval,
            s1 * dq1 * (w_bar_eval - 1) - s2 * dq2 * w_bar_eval,
            -cos_dq * ddq_diff * h / D,
        ])

        # en = en_unnorm / delta => den/dt = (1/delta)*(den_unnorm/dt)
        #                                   - (en_unnorm/delta^2)*(ddelta/dt)
        d_delta_sq = (
            2 * D**2 * w_bar_eval * (w_bar_eval - 1) * sin_dq * ddq_diff
            + 2 * h**2 * sin_dq * cos_dq * ddq_diff
        )
        d_delta = d_delta_sq / (2 * delta) if abs(delta) > 1e-30 else 0.0

        en_unnorm_val = D * np.array([
            s1 * (w_bar_eval - 1) - s2 * w_bar_eval,
            -c1 * (w_bar_eval - 1) + c2 * w_bar_eval,
            -sin_dq * h / D,
        ])

        den_dt = den_unnorm / delta - en_unnorm_val * d_delta / delta**2
        return den_dt

    def _vn_dot_at(h, w_bar_eval):
        return _compute_vn_dot(
            h, w_bar_eval, q1, q2, dq1, dq2, ddq1, ddq2, D,
        )

    # --- Trailing-edge shedding term: 0.5 * m(hmax) * vn^2 * en ---
    # Evaluate at the midpoint w_bar=0.5 (representative strip center)
    w_bar_mid = 0.5
    vn_te, en_te = _vn_en_at(hmax, w_bar_mid)
    m_te = _added_mass_per_h(hmax)
    F_shed = 0.5 * m_te * vn_te**2 * en_te

    # --- Bulk momentum rate: -d/dt integral m(h)*vn*en dh ---
    # Full product rule: d/dt[m*vn*en] = (dm/dt)*vn*en + m*(dvn/dt)*en + m*vn*(den/dt)
    from scipy.integrate import quad

    F_momentum_rate = np.zeros(3)
    for k in range(3):
        def _momentum_integrand(h, _k=k):
            m_h = _added_mass_per_h(h)
            dm_h = _dm_dt_per_h(h)
            vn, en = _vn_en_at(h, w_bar_mid)
            vn_dot = _vn_dot_at(h, w_bar_mid)
            den_dt = _den_dt_at(h, w_bar_mid)
            # Full product rule: dm/dt*vn*en + m*dvn/dt*en + m*vn*den/dt
            return (dm_h * vn * en[_k]
                    + m_h * vn_dot * en[_k]
                    + m_h * vn * den_dt[_k])

        val, _ = quad(_momentum_integrand, hmin, hmax, epsabs=1e-8, epsrel=1e-8)
        F_momentum_rate[k] = val

    # Lighthill's reactive force = -momentum_rate + shedding
    F_reactive = -F_momentum_rate + F_shed

    return F_reactive


def compute_panel_torques(
    D: float,
    hmin: float,
    hmax: float,
    q1: float,
    q2: float,
    dq1: float,
    dq2: float,
    Cn: float,
    rho: float,
) -> tuple[float, float]:
    """Compute hydrodynamic torques on joints from a membrane panel.

    Thesis Eq. (15) context — torques τ_{Φ,i} on joints due to hydrodynamic
    loads. The torque on each joint from panel Ω_{1|2} is obtained by
    integrating the moment arm × normal force over the membrane.

    The torque on joint 1 uses the lever arm h·cos(q1)·(w̄−1) and
    the torque on joint 2 uses h·cos(q2)·w̄/w_max.

    Parameters
    ----------
    D, hmin, hmax, q1, q2, dq1, dq2, Cn, rho : float
        Same as compute_panel_forces.

    Returns
    -------
    tau1 : float
        Torque on the first joint of the panel [N·m].
    tau2 : float
        Torque on the second joint of the panel [N·m].
    """
    dq = q1 - q2
    cos_dq = np.cos(dq)
    sin_dq = np.sin(dq)
    s1, c1 = np.sin(q1), np.cos(q1)
    s2, c2 = np.sin(q2), np.cos(q2)

    def w_max(h: float) -> float:
        return np.sqrt(2 * h**2 * (1 - cos_dq) + D**2)

    def _torque_integrand(w: float, h: float, joint: int) -> float:
        wm = w_max(h)
        if wm < 1e-30:
            return 0.0
        w_bar = w / wm

        delta = np.sqrt(
            2 * D**2 * w_bar * (w_bar - 1) * (1 - cos_dq)
            + h**2 * sin_dq**2
            + D**2
        )
        if abs(delta) < 1e-30:
            return 0.0

        vn = (D * h / delta) * (
            dq1
            + (1 - cos_dq) * (dq1 + dq2) * w_bar**2
            + (cos_dq * (dq1 + dq2) - 2 * dq1) * w_bar
        )

        sgn_vn = np.sign(vn) if abs(vn) > 1e-30 else 0.0

        # For torque, we need the dot product of (moment_arm × f_n)
        # projected along the joint axis. The heave (en_y) term contributes
        # to joint torque through cos(q_i) lever and sway (en_x) through sin(q_i).
        en_x = (D / delta) * (s1 * (w_bar - 1) - s2 * w_bar)
        en_y = (D / delta) * (-c1 * (w_bar - 1) + c2 * w_bar)

        f_n_mag = -0.5 * rho * Cn * vn**2 * sgn_vn

        if joint == 0:
            # Torque on joint 1: lever arm = h·(w̄−1) factor
            lever_x = c1 * (w_bar - 1)
            lever_y = s1 * (w_bar - 1)
            torque_contrib = h * D * (lever_x * en_y - lever_y * en_x) * f_n_mag
        else:
            # Torque on joint 2: lever arm with w̄ and 1/w_max
            lever_x = -c2 * w_bar / wm
            lever_y = -s2 * w_bar / wm
            torque_contrib = h * D * wm * (lever_x * en_y - lever_y * en_x) * f_n_mag

        return torque_contrib

    tau1, _ = dblquad(
        lambda w, h: _torque_integrand(w, h, 0),
        hmin, hmax, 0, w_max,
        epsabs=1e-8, epsrel=1e-8,
    )
    tau2, _ = dblquad(
        lambda w, h: _torque_integrand(w, h, 1),
        hmin, hmax, 0, w_max,
        epsabs=1e-8, epsrel=1e-8,
    )

    return tau1, tau2
