"""The prescribed motion must carry the reference surface onto the
surface at the new time, exactly.

If it does not, the CFD is simulating a differently-shaped body at
every instant and no comparison with the analytical model is valid.
"""

import numpy as np
import pytest

from undulating_fin.params import FinParams
from openfoam.geometry import (
    SPAN_SIGN,
    membrane_displacement,
    ray_angles,
    surface_point,
)


PARAMS = [
    FinParams(N=2, D=0.09, hmin=0.034, hmax=0.174, A=np.deg2rad(30),
              f=1.0, Phi=np.deg2rad(0)),
    FinParams(N=3, D=0.09, hmin=0.034, hmax=0.174, A=np.deg2rad(30),
              f=1.0, Phi=np.deg2rad(40)),
    FinParams(N=6, D=0.10, hmin=0.0, hmax=0.175, A=np.deg2rad(25),
              f=2.0, Phi=np.deg2rad(60), X=0.05),
]
TIMES = [0.0, 0.13, 0.25, 0.5, 0.77, 1.0]


def _sample_surface(p: FinParams, t: float, n_h=7, n_w=5):
    """Mid-surface sample points and their (panel, w̄, h) labels."""
    q = ray_angles(p, t)
    pts, labels = [], []
    for panel in range(p.N - 1):
        z0 = SPAN_SIGN * panel * p.D
        for h in np.linspace(p.hmin, p.hmax, n_h):
            for w in np.linspace(0.0, 1.0, n_w):
                pts.append(surface_point(h, w, q[panel], q[panel + 1], p.D, z0))
                labels.append((panel, w, h))
    return np.asarray(pts), labels


@pytest.mark.parametrize("p", PARAMS, ids=lambda p: f"N{p.N}")
@pytest.mark.parametrize("t", TIMES)
def test_displacement_maps_reference_onto_current_surface(p, t):
    """P₀ + d(t) must equal the analytic surface point at time t."""
    p0, labels = _sample_surface(p, 0.0)
    moved = p0 + membrane_displacement(p0, p, t, t_ref=0.0)

    q = ray_angles(p, t)
    expected = np.array([
        surface_point(h, w, q[panel], q[panel + 1], p.D,
                      SPAN_SIGN * panel * p.D)
        for (panel, w, h) in labels
    ])
    assert moved == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("p", PARAMS, ids=lambda p: f"N{p.N}")
def test_zero_displacement_at_reference_time(p):
    p0, _ = _sample_surface(p, 0.0)
    d = membrane_displacement(p0, p, t=0.0, t_ref=0.0)
    assert np.abs(d).max() < 1e-12


@pytest.mark.parametrize("p", PARAMS, ids=lambda p: f"N{p.N}")
@pytest.mark.parametrize("t", TIMES)
def test_motion_is_purely_in_plane(p, t):
    """Rays rotate about z, so the span coordinate never changes."""
    p0, _ = _sample_surface(p, 0.0)
    assert np.abs(membrane_displacement(p0, p, t)[..., 2]).max() == 0.0


@pytest.mark.parametrize("p", PARAMS, ids=lambda p: f"N{p.N}")
@pytest.mark.parametrize("t", TIMES)
def test_points_on_a_ray_keep_their_radius(p, t):
    """Points *on* a ray are rigidly rotated, so their radius is fixed.

    Interior membrane points are not: they sit on the straight chord
    between two rays, and the chord shortens as the rays splay apart.
    Only w̄ = 0 and w̄ = 1 lie on a ray, where g = 1 and r = h.
    """
    q0 = ray_angles(p, 0.0)
    q = ray_angles(p, t)
    for panel in range(p.N - 1):
        z0 = SPAN_SIGN * panel * p.D
        for w in (0.0, 1.0):
            for h in np.linspace(p.hmin, p.hmax, 5):
                pt = surface_point(h, w, q0[panel], q0[panel + 1], p.D, z0)
                moved = pt + membrane_displacement(pt, p, t)
                assert np.hypot(moved[0], moved[1]) == pytest.approx(h, rel=1e-12)
                assert np.hypot(pt[0], pt[1]) == pytest.approx(h, rel=1e-12)


@pytest.mark.parametrize("p", PARAMS, ids=lambda p: f"N{p.N}")
@pytest.mark.parametrize("t", TIMES)
def test_interior_radius_follows_the_chord(p, t):
    """Off a ray the radius must track r = h·g(w̄, Δq(t)) exactly.

    Pins the one non-obvious consequence of the ruled geometry: the
    membrane's mid-chord is pulled *inwards* as the bounding rays
    splay, so a constant-radius assumption is wrong.
    """
    q0 = ray_angles(p, 0.0)
    q = ray_angles(p, t)
    for panel in range(p.N - 1):
        z0 = SPAN_SIGN * panel * p.D
        dq = q[panel] - q[panel + 1]
        for w in (0.25, 0.5, 0.75):
            g = np.sqrt(1.0 + 2.0 * w * (w - 1.0) * (1.0 - np.cos(dq)))
            for h in np.linspace(max(p.hmin, 1e-3), p.hmax, 4):
                pt = surface_point(h, w, q0[panel], q0[panel + 1], p.D, z0)
                moved = pt + membrane_displacement(pt, p, t)
                assert np.hypot(moved[0], moved[1]) == pytest.approx(
                    h * g, rel=1e-12
                )


@pytest.mark.parametrize("p", PARAMS, ids=lambda p: f"N{p.N}")
def test_periodicity(p):
    """Motion must repeat exactly after one oscillation period."""
    p0, _ = _sample_surface(p, 0.0)
    T = 1.0 / p.f
    for t in (0.2, 0.6):
        assert membrane_displacement(p0, p, t) == pytest.approx(
            membrane_displacement(p0, p, t + T), abs=1e-12
        )


@pytest.mark.parametrize("t", TIMES)
def test_inflated_solid_survives_the_motion(t):
    """The thickened STL must move without collapsing or inverting.

    Most mesh points are off the mid-surface, so the map has to carry
    the whole solid, not just the sheet.  Checks the solid keeps a
    positive volume of the right order and that no two vertices are
    driven onto each other.
    """
    from openfoam.stl import build_fin_solid, signed_volume

    p = FinParams(N=3, D=0.09, hmin=0.034, hmax=0.174,
                  A=np.deg2rad(30), f=1.0, Phi=np.deg2rad(40))
    verts, faces = build_fin_solid(p, t=0.0, n_h=12, n_w=6, thickness=4e-3)
    v0 = signed_volume(verts, faces)

    moved = verts + membrane_displacement(verts, p, t)
    vt = signed_volume(moved, faces)

    assert vt > 0.0
    assert vt == pytest.approx(v0, rel=0.10)
    assert moved[:, 2] == pytest.approx(verts[:, 2], abs=1e-15)

    # No vertex collapse: the closest distinct pair stays well separated.
    from scipy.spatial import cKDTree
    d, _ = cKDTree(moved).query(moved, k=2)
    assert d[:, 1].min() > 1e-4
