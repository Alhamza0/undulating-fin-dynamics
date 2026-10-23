r"""
Membrane surface geometry for the undulating fin.

This module defines the *exact* parametric surface that the analytical
model integrates over, so that the CFD mesh and the analytical model
describe the same body.  Everything here is derived from the thesis
equations; :func:`surface_normal_and_delta` is verified numerically
against :func:`undulating_fin.hydrodynamics._compute_vn_en_delta`
(Eq. 15 and Eq. 16) by the test-suite.

Parameterisation
----------------
A membrane panel spans two adjacent rays.  Ray *i* lies at spanwise
station ``z = −i·D`` and is rotated by its joint angle ``q_i(t)`` about
the *z* axis; the membrane between two rays is the ruled surface

.. math::

    \mathbf{P}(h, \bar w) = \big(\,
        h[(1-\bar w)\cos q_1 + \bar w \cos q_2],\;
        h[(1-\bar w)\sin q_1 + \bar w \sin q_2],\;
        -\bar w D \,\big)

with ``h ∈ [hmin, hmax]`` the height along the ray and
``w̄ ∈ [0, 1]`` the normalised spanwise coordinate.

.. warning::

    The rays stack along **−z**, not +z.  This is not cosmetic.  The
    surface normal works out to

    .. math::

        \partial_h \mathbf{P} \times \partial_{\bar w} \mathbf{P}
        = \big(s D V,\; -s D U,\; -h \sin \Delta q\big)

    with :math:`U = (1-\bar w)\cos q_1 + \bar w \cos q_2` and
    :math:`V = (1-\bar w)\sin q_1 + \bar w \sin q_2`, and *s* is
    ``SPAN_SIGN``.  Reversing the span therefore flips the
    *x* and *y* components while leaving *z* untouched — i.e. it
    silently inverts **heave and sway** but not surge.  δ, and hence the
    surface area, is unchanged either way, so nothing else gives the
    error away.  ``SPAN_SIGN`` below fixes the convention and
    ``test_normal_matches_thesis_eq16`` pins it to Eq. (16).

Taking the cross product of the tangents reproduces Eq. (15) exactly:

.. math::

    \delta = \left|
        \partial_h \mathbf{P} \times \partial_{\bar w} \mathbf{P}
    \right|
    = \sqrt{2D^2\bar w(\bar w-1)(1-\cos\Delta q)
             + h^2\sin^2\!\Delta q + D^2}

which is the geometric scale factor used by the force integral.
"""

from __future__ import annotations

import numpy as np

from undulating_fin.params import FinParams

#: Spanwise stacking direction of the rays.  See the module warning:
#: this sign controls the heave/sway force orientation and must stay −1
#: to match thesis Eq. (16).
SPAN_SIGN = -1.0

__all__ = [
    "SPAN_SIGN",
    "ray_angles",
    "surface_point",
    "surface_normal_and_delta",
    "membrane_mesh",
    "membrane_displacement",
]


# ── Kinematics ──────────────────────────────────────────────────────

def ray_angles(params: FinParams, t: float) -> np.ndarray:
    """Joint angles ``q_i(t)`` for all rays [rad].  Thesis Eq. (1).

    ``q_i(t) = A·sin(2πft − i·Φ) + X`` with *i* zero-indexed.
    """
    i = np.arange(params.N)
    return params.A * np.sin(2 * np.pi * params.f * t - i * params.Phi) + params.X


# ── Surface ─────────────────────────────────────────────────────────

def surface_point(
    h: np.ndarray | float,
    w_bar: np.ndarray | float,
    q1: float,
    q2: float,
    D: float,
    z0: float = 0.0,
) -> np.ndarray:
    """Point on the ruled membrane surface.

    Parameters
    ----------
    h
        Height along the ray [m].
    w_bar
        Normalised spanwise coordinate in ``[0, 1]``.
    q1, q2
        Joint angles of the bounding rays [rad].
    D
        Inter-ray distance [m].
    z0
        Spanwise station of the first ray of this panel [m], normally
        ``SPAN_SIGN · panel_index · D``.

    Returns
    -------
    ndarray
        Shape ``(..., 3)`` array of ``(x, y, z)`` coordinates.
    """
    h = np.asarray(h, dtype=float)
    w_bar = np.asarray(w_bar, dtype=float)

    # Radial unit direction, linearly blended between the two rays.
    ux = (1.0 - w_bar) * np.cos(q1) + w_bar * np.cos(q2)
    uy = (1.0 - w_bar) * np.sin(q1) + w_bar * np.sin(q2)

    return np.stack([h * ux, h * uy, z0 + SPAN_SIGN * w_bar * D], axis=-1)


def surface_normal_and_delta(
    h: float,
    w_bar: float,
    q1: float,
    q2: float,
    D: float,
) -> tuple[np.ndarray, float]:
    """Unit normal and scale factor ``δ`` from the parameterisation.

    Computed as ``∂P/∂h × ∂P/∂w̄`` — i.e. from first principles rather
    than from the closed-form thesis expression.  Agreement between the
    two is what the test-suite checks.

    Returns
    -------
    en : ndarray, shape (3,)
        Unit normal vector.
    delta : float
        ``|∂P/∂h × ∂P/∂w̄|`` — thesis Eq. (15) [m].
    """
    c1, c2 = np.cos(q1), np.cos(q2)
    s1, s2 = np.sin(q1), np.sin(q2)

    # ∂P/∂h — the blended radial direction (independent of h).
    dP_dh = np.array([
        (1.0 - w_bar) * c1 + w_bar * c2,
        (1.0 - w_bar) * s1 + w_bar * s2,
        0.0,
    ])

    # ∂P/∂w̄
    dP_dw = np.array([
        h * (c2 - c1),
        h * (s2 - s1),
        SPAN_SIGN * D,
    ])

    n = np.cross(dP_dh, dP_dw)
    delta = float(np.linalg.norm(n))
    if delta < 1e-30:
        return np.zeros(3), 1e-30
    return n / delta, delta


# ── Discretisation ──────────────────────────────────────────────────

def membrane_mesh(
    params: FinParams,
    t: float = 0.0,
    n_h: int = 40,
    n_w: int = 20,
) -> tuple[np.ndarray, np.ndarray]:
    """Tessellate the whole fin (all panels) into a triangle mesh.

    Parameters
    ----------
    params
        Fin parameters.
    t
        Time at which to evaluate the ray angles [s].
    n_h, n_w
        Number of subdivisions along the height and along each panel's
        span.  Vertex counts are ``n_h + 1`` and ``n_w + 1``.

    Returns
    -------
    vertices : ndarray, shape (n_vertices, 3)
    faces : ndarray, shape (n_faces, 3)
        Triangle vertex indices, wound counter-clockwise as seen from
        the ``+normal`` side.
    """
    q = ray_angles(params, t)
    h_vals = np.linspace(params.hmin, params.hmax, n_h + 1)
    w_vals = np.linspace(0.0, 1.0, n_w + 1)

    verts: list[np.ndarray] = []
    faces: list[tuple[int, int, int]] = []

    for panel in range(params.N - 1):
        base = len(verts)
        q1, q2 = q[panel], q[panel + 1]
        z0 = SPAN_SIGN * panel * params.D

        # Grid of points, indexed [ih, iw].
        for ih, h in enumerate(h_vals):
            for iw, w in enumerate(w_vals):
                verts.append(surface_point(h, w, q1, q2, params.D, z0))

        stride = n_w + 1
        for ih in range(n_h):
            for iw in range(n_w):
                v00 = base + ih * stride + iw
                v01 = v00 + 1
                v10 = v00 + stride
                v11 = v10 + 1
                faces.append((v00, v01, v11))
                faces.append((v00, v11, v10))

    return np.asarray(verts, dtype=float), np.asarray(faces, dtype=int)


# ── Prescribed motion ───────────────────────────────────────────────

def membrane_displacement(
    p0: np.ndarray,
    params: FinParams,
    t: float,
    t_ref: float = 0.0,
) -> np.ndarray:
    r"""Displacement of membrane points from their reference positions.

    This is the reference implementation of the motion imposed by the
    ``pointDisplacement`` boundary condition; the C++ in the generated
    case is a transcription of it, and
    ``test_cpp_transcription_matches_python`` keeps the two in step.

    A mesh point is located by its spanwise station, since
    ``z = -(panel + w̄)·D`` gives the panel and w̄ directly:

    .. math::
        s = -z_0 / D, \quad
        \text{panel} = \lfloor s \rfloor, \quad
        \bar w = s - \text{panel}

    The blended radial direction
    :math:`\mathbf u = (1-\bar w)\hat u(q_1) + \bar w \hat u(q_2)`
    is *not* a unit vector once the bounding rays differ; it has length
    :math:`g = \sqrt{1 + 2\bar w(\bar w - 1)(1 - \cos \Delta q)}`
    and direction :math:`\varphi = \operatorname{atan2}(u_y, u_x)`.
    The motion of the station is therefore a rotation by
    :math:`\Delta\varphi` together with a radial stretch by
    :math:`g(t)/g(t_{\mathrm{ref}})`:

    .. math::
        \mathbf x_1 = \frac{g(t)}{g(t_{\mathrm{ref}})}
                       R(\Delta\varphi)\, \mathbf x_0

    Both factors depend only on the station, never on where the point
    sits within it, which is what lets the map carry the *inflated*
    solid rather than just the sheet.  Recovering ``h`` from the radius
    and re-evaluating the surface would instead project every point
    onto the mid-surface and collapse the membrane's thickness.

    Rays only rotate about *z*, so the span coordinate is untouched and
    the displacement lies wholly in the *x*-*y* plane.

    Parameters
    ----------
    p0
        Reference positions, shape ``(..., 3)``.
    params
        Fin parameters.
    t
        Current time [s].
    t_ref
        Time at which the reference configuration was built [s].

    Returns
    -------
    ndarray
        Displacements, same shape as *p0*.
    """
    p0 = np.asarray(p0, dtype=float)
    x0, y0, z0 = p0[..., 0], p0[..., 1], p0[..., 2]

    s = -z0 / params.D
    panel = np.clip(np.floor(s).astype(int), 0, max(params.N - 2, 0))
    w = np.clip(s - panel, 0.0, 1.0)

    q_ref = ray_angles(params, t_ref)
    q_now = ray_angles(params, t)

    def blend(q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return (
            (1.0 - w) * np.cos(q[panel]) + w * np.cos(q[panel + 1]),
            (1.0 - w) * np.sin(q[panel]) + w * np.sin(q[panel + 1]),
        )

    ux_r, uy_r = blend(q_ref)
    ux_t, uy_t = blend(q_now)

    g_ref = np.maximum(np.hypot(ux_r, uy_r), 1e-300)
    g_now = np.hypot(ux_t, uy_t)

    # Rotate by the swing of the blended direction and stretch by the
    # change in its length.  Both are properties of the station (w̄, t)
    # alone, so the map is well defined off the mid-surface too.
    dphi = np.arctan2(uy_t, ux_t) - np.arctan2(uy_r, ux_r)
    scale = g_now / g_ref
    cd, sd = np.cos(dphi), np.sin(dphi)

    x1 = scale * (x0 * cd - y0 * sd)
    y1 = scale * (x0 * sd + y0 * cd)

    return np.stack([x1 - x0, y1 - y0, np.zeros_like(z0)], axis=-1)
