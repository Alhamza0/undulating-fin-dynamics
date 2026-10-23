"""Verify the OpenFOAM surface parameterisation against the thesis equations.

The CFD mesh is built from :mod:`openfoam.geometry`.  The analytical
force model integrates over the surface described by Eq. (15)/(16) in
:mod:`undulating_fin.hydrodynamics`.  If the two disagree, the CFD and
the analytical model are simulating different bodies and any comparison
between them is meaningless.  These tests pin them together.
"""

import numpy as np
import pytest

from undulating_fin.params import FinParams
from undulating_fin.hydrodynamics import _compute_vn_en_delta
from openfoam.geometry import (
    SPAN_SIGN,
    ray_angles,
    surface_point,
    surface_normal_and_delta,
    membrane_mesh,
)


# Representative sample of the (h, w̄, q1, q2) space, including the
# degenerate q1 == q2 case where the panel is planar.
CASES = [
    (0.05, 0.0, 0.0, 0.0),
    (0.10, 0.5, 0.3, -0.2),
    (0.174, 1.0, 0.5236, -0.5236),
    (0.034, 0.25, -0.4, 0.7),
    (0.12, 0.75, 0.1, 0.1),
    (0.15, 0.33, 1.0, -0.9),
]
D = 0.09


@pytest.mark.parametrize("h,w,q1,q2", CASES)
def test_delta_matches_thesis_eq15(h, w, q1, q2):
    """|∂P/∂h × ∂P/∂w̄| must equal the closed-form δ of Eq. (15)."""
    _, _, delta_thesis = _compute_vn_en_delta(h, w, q1, q2, 0.0, 0.0, D)
    _, delta_geom = surface_normal_and_delta(h, w, q1, q2, D)
    assert delta_geom == pytest.approx(delta_thesis, rel=1e-12)


@pytest.mark.parametrize("h,w,q1,q2", CASES)
def test_normal_matches_thesis_eq16(h, w, q1, q2):
    """The parameterisation's unit normal must equal Eq. (16) exactly.

    Not just collinear: the *signed* match is what guarantees the CFD
    surge force has the same orientation as the analytical one.  This
    only holds with ``SPAN_SIGN == -1`` (rays stacking along −z); with
    the span reversed the z component flips and thrust changes sign.
    """
    _, en_thesis, _ = _compute_vn_en_delta(h, w, q1, q2, 0.0, 0.0, D)
    en_geom, _ = surface_normal_and_delta(h, w, q1, q2, D)
    assert np.linalg.norm(en_thesis) == pytest.approx(1.0, rel=1e-12)
    assert en_geom == pytest.approx(en_thesis, abs=1e-12)


@pytest.mark.parametrize("h,w,q1,q2", CASES)
def test_reversed_span_inverts_heave_and_sway_only(h, w, q1, q2):
    """Guard rail documenting exactly what a wrong span direction costs.

    Reversing SPAN_SIGN flips the normal's x (heave) and y (sway)
    components and leaves z (surge) alone.  Building the CFD mesh with
    the span the wrong way round would therefore invert two of the
    three reported force components while δ — and so the fin's surface
    area — looks perfectly correct.
    """
    import openfoam.geometry as geo

    _, en_thesis, _ = _compute_vn_en_delta(h, w, q1, q2, 0.0, 0.0, D)
    original = geo.SPAN_SIGN
    try:
        geo.SPAN_SIGN = -original
        en_flipped, _ = geo.surface_normal_and_delta(h, w, q1, q2, D)
    finally:
        geo.SPAN_SIGN = original

    assert en_flipped[0] == pytest.approx(-en_thesis[0], abs=1e-12)
    assert en_flipped[1] == pytest.approx(-en_thesis[1], abs=1e-12)
    assert en_flipped[2] == pytest.approx(+en_thesis[2], abs=1e-12)


@pytest.mark.parametrize("h,w,q1,q2", CASES)
def test_delta_is_the_jacobian_of_the_surface(h, w, q1, q2):
    """δ must equal the finite-difference surface Jacobian.

    Independent check that δ really is the area scale factor dS/(dh dw̄),
    not just an algebraic coincidence.
    """
    eps = 1e-7
    p0 = surface_point(h, w, q1, q2, D)
    dh = (surface_point(h + eps, w, q1, q2, D) - p0) / eps
    dw = (surface_point(h, w + eps, q1, q2, D) - p0) / eps
    assert np.linalg.norm(np.cross(dh, dw)) == pytest.approx(
        surface_normal_and_delta(h, w, q1, q2, D)[1], rel=1e-5
    )


def test_ray_angles_match_eq1():
    """Eq. (1): q_i(t) = A·sin(2πft − iΦ) + X."""
    p = FinParams(N=4, A=np.deg2rad(30), f=1.0, Phi=np.deg2rad(40), X=0.1)
    t = 0.37
    expected = [
        p.A * np.sin(2 * np.pi * p.f * t - i * p.Phi) + p.X for i in range(p.N)
    ]
    assert ray_angles(p, t) == pytest.approx(expected, rel=1e-14)


def test_rays_sit_at_correct_spanwise_stations():
    """w̄ = 0 and w̄ = 1 must land exactly on the two bounding rays."""
    p = FinParams(N=3, D=0.09, hmin=0.034, hmax=0.174)
    q = ray_angles(p, 0.0)
    for panel in range(p.N - 1):
        z0 = SPAN_SIGN * panel * p.D
        a = surface_point(p.hmax, 0.0, q[panel], q[panel + 1], p.D, z0)
        b = surface_point(p.hmax, 1.0, q[panel], q[panel + 1], p.D, z0)
        # Spanwise stations
        assert a[2] == pytest.approx(z0)
        assert b[2] == pytest.approx(z0 + SPAN_SIGN * p.D)
        # Radial distance equals h exactly on the rays
        assert np.hypot(a[0], a[1]) == pytest.approx(p.hmax)
        assert np.hypot(b[0], b[1]) == pytest.approx(p.hmax)
        # And the ray makes angle q with the x axis
        assert np.arctan2(a[1], a[0]) == pytest.approx(q[panel])
        assert np.arctan2(b[1], b[0]) == pytest.approx(q[panel + 1])


def test_membrane_mesh_is_conformal_and_finite():
    p = FinParams(N=3, D=0.09, hmin=0.034, hmax=0.174)
    verts, faces = membrane_mesh(p, t=0.0, n_h=10, n_w=5)
    n_panels = p.N - 1
    assert len(verts) == n_panels * 11 * 6
    assert len(faces) == n_panels * 10 * 5 * 2
    assert np.all(np.isfinite(verts))
    assert faces.max() == len(verts) - 1
    # No degenerate triangles
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    assert areas.min() > 0.0


def test_total_area_converges_to_the_force_integral_domain():
    """Summed triangle area must converge to ∬ δ dh dw̄ — the exact
    surface area the analytical force integral runs over."""
    from scipy.integrate import dblquad

    p = FinParams(N=2, D=0.09, hmin=0.034, hmax=0.174)
    q = ray_angles(p, 0.31)

    exact, _ = dblquad(
        lambda w, h: surface_normal_and_delta(h, w, q[0], q[1], p.D)[1],
        p.hmin, p.hmax, 0.0, 1.0, epsabs=1e-12,
    )

    verts, faces = membrane_mesh(p, t=0.31, n_h=200, n_w=100)
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    tessellated = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()

    assert tessellated == pytest.approx(exact, rel=1e-4)
