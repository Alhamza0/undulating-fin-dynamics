"""The fin STL must be a closed, correctly-oriented solid.

snappyHexMesh decides inside from outside by ray casting, so a leaky or
inside-out surface produces a mesh that looks plausible and is wrong.
"""

import numpy as np
import pytest

from undulating_fin.params import FinParams
from openfoam.stl import (
    build_fin_solid,
    is_watertight,
    orient_outward,
    signed_volume,
    write_stl_ascii,
    write_stl_binary,
)
from openfoam.geometry import ray_angles, surface_point


@pytest.fixture(params=[2, 3, 5])
def solid(request):
    p = FinParams(N=request.param, D=0.09, hmin=0.034, hmax=0.174)
    return p, *build_fin_solid(p, t=0.0, n_h=12, n_w=6, thickness=2e-3)


def test_watertight(solid):
    _, _, faces = solid
    assert is_watertight(faces)


def test_no_degenerate_triangles(solid):
    _, verts, faces = solid
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    assert areas.min() > 1e-12


def test_closed_surface_encloses_positive_volume(solid):
    """Divergence theorem: V = (1/6)Σ a·(b×c).

    Positive volume proves the winding is outward everywhere — an
    inside-out solid gives exactly the negative of the right answer.
    """
    p, verts, faces = solid
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    volume = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0
    assert volume > 0.0

    # Should be close to (membrane area) × thickness.
    span = (p.N - 1) * p.D
    approx_area = span * (p.hmax - p.hmin)
    assert volume == pytest.approx(approx_area * 2e-3, rel=0.15)


def test_zero_net_area_vector(solid):
    """Σ (facet area × normal) must vanish for any closed surface."""
    _, verts, faces = solid
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    assert np.allclose(0.5 * np.cross(b - a, c - a).sum(axis=0), 0.0, atol=1e-12)


def test_hull_brackets_the_mid_surface(solid):
    """Every mid-surface point must lie inside the inflated solid's
    bounding box — i.e. the offset really did wrap the sheet."""
    p, verts, _ = solid
    q = ray_angles(p, 0.0)
    lo, hi = verts.min(axis=0), verts.max(axis=0)
    for panel in range(p.N - 1):
        for w in (0.0, 0.5, 1.0):
            for h in (p.hmin, p.hmax):
                pt = surface_point(h, w, q[panel], q[panel + 1], p.D,
                                   -panel * p.D)
                assert np.all(pt >= lo - 1e-12) and np.all(pt <= hi + 1e-12)


def test_solid_occupies_the_expected_envelope(solid):
    """Bounding box must match the fin's actual reach.

    Radially the sheet cannot exceed hmax; spanwise it must cover
    exactly (N-1)·D, running along -z per SPAN_SIGN.
    """
    p, verts, _ = solid
    span = (p.N - 1) * p.D
    radial = np.hypot(verts[:, 0], verts[:, 1])

    assert radial.max() <= p.hmax + 2e-3
    assert radial.min() >= p.hmin - 2e-3
    assert verts[:, 2].max() == pytest.approx(0.0, abs=2e-3)
    assert verts[:, 2].min() == pytest.approx(-span, abs=2e-3)


def test_invalid_thickness_rejected():
    p = FinParams(N=2)
    with pytest.raises(ValueError):
        build_fin_solid(p, thickness=0.0)


def test_writers_roundtrip(tmp_path, solid):
    _, verts, faces = solid

    a = write_stl_ascii(tmp_path / "fin.stl", verts, faces)
    txt = a.read_text()
    assert txt.startswith("solid fin")
    assert txt.count("facet normal") == len(faces)

    b = write_stl_binary(tmp_path / "fin_bin.stl", verts, faces)
    raw = b.read_bytes()
    assert len(raw) == 84 + 50 * len(faces)
    assert int.from_bytes(raw[80:84], "little") == len(faces)


def test_orientation_pass_repairs_a_scrambled_winding(solid):
    """orient_outward must recover from arbitrarily flipped facets.

    This is the property the STL builder relies on instead of a
    hand-derived winding for every shell and rim.
    """
    _, verts, faces = solid
    rng = np.random.default_rng(0)
    scrambled = faces.copy()
    flip = rng.random(len(faces)) < 0.5
    scrambled[flip] = scrambled[flip][:, ::-1]

    repaired = orient_outward(verts, scrambled)

    assert signed_volume(verts, repaired) == pytest.approx(
        signed_volume(verts, faces), rel=1e-12
    )
    a, b, c = verts[repaired[:, 0]], verts[repaired[:, 1]], verts[repaired[:, 2]]
    assert np.allclose(0.5 * np.cross(b - a, c - a).sum(axis=0), 0.0, atol=1e-12)


def test_orientation_pass_fixes_a_fully_inverted_solid(solid):
    """An inside-out solid must come back with positive volume."""
    _, verts, faces = solid
    inverted = faces[:, ::-1]
    assert signed_volume(verts, inverted) < 0.0
    assert signed_volume(verts, orient_outward(verts, inverted)) > 0.0
