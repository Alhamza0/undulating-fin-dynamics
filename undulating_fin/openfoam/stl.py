r"""
Watertight STL generation for the undulating fin membrane.

``snappyHexMesh`` needs a closed surface, but the membrane of thesis
Eq. (15) is a zero-thickness ruled sheet.  This module inflates it into
a thin closed solid: the mid-surface is offset by ``±thickness/2`` along
its unit normal and the two shells are stitched together around the
boundary.

The whole fin is built as **one** grid, so panels that meet at a ray
share their vertices exactly and the result stays watertight for any
number of rays.

Winding follows the STL convention — counter-clockwise seen from
outside, so the facet normal by the right-hand rule points out of the
solid.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

from undulating_fin.params import FinParams

from .geometry import SPAN_SIGN, ray_angles, surface_normal_and_delta, surface_point

__all__ = [
    "build_fin_solid",
    "write_stl_ascii",
    "write_stl_binary",
    "is_watertight",
    "orient_outward",
    "signed_volume",
]


# ── Orientation ─────────────────────────────────────────────────────

def signed_volume(verts: np.ndarray, faces: np.ndarray) -> float:
    """Signed volume of a closed surface, ``(1/6)Σ a·(b×c)``.

    Positive when the winding is outward.  Meaningless for an open or
    inconsistently-wound surface.
    """
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def orient_outward(verts: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Return *faces* rewound consistently, with normals pointing out.

    Hand-deriving the winding of every shell and rim is a reliable
    source of sign slips, and a single inverted rim is invisible in a
    rendering but makes ``snappyHexMesh`` carve the mesh inside out.
    So the winding is repaired here instead:

    1. Two triangles sharing an edge agree iff they traverse that edge
       in *opposite* directions.  Flood-filling from one seed face and
       flipping on disagreement makes the whole surface consistent.
    2. That leaves a global two-fold ambiguity, settled by the sign of
       the enclosed volume.

    Requires a watertight surface — check with :func:`is_watertight`.
    """
    faces = np.array(faces, dtype=int, copy=True)

    # Undirected edge -> the faces touching it.
    edge_faces: dict[tuple[int, int], list[int]] = {}
    for fi, (a, b, c) in enumerate(faces):
        for u, v in ((a, b), (b, c), (c, a)):
            edge_faces.setdefault((u, v) if u < v else (v, u), []).append(fi)

    def directed(fi: int) -> set[tuple[int, int]]:
        a, b, c = faces[fi]
        return {(a, b), (b, c), (c, a)}

    seen = np.zeros(len(faces), dtype=bool)
    seen[0] = True
    stack = [0]
    while stack:
        fi = stack.pop()
        for e in tuple(directed(fi)):
            key = e if e[0] < e[1] else (e[1], e[0])
            for fj in edge_faces.get(key, ()):
                if fj == fi or seen[fj]:
                    continue
                # Same direction on the shared edge => fj is flipped.
                if e in directed(fj):
                    faces[fj] = faces[fj][::-1]
                seen[fj] = True
                stack.append(fj)

    if signed_volume(verts, faces) < 0.0:
        faces = faces[:, ::-1]
    return faces


# ── Surface construction ────────────────────────────────────────────

def _mid_surface(
    params: FinParams, t: float, n_h: int, n_w: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Mid-surface grid and its unit normals.

    Returns arrays of shape ``(n_h + 1, n_span + 1, 3)`` where
    ``n_span = (N - 1) * n_w``.  Grid columns that fall on a ray are
    shared between the adjacent panels, so the sheet is continuous.
    """
    q = ray_angles(params, t)
    h_vals = np.linspace(params.hmin, params.hmax, n_h + 1)
    n_span = (params.N - 1) * n_w

    pts = np.empty((n_h + 1, n_span + 1, 3))
    nrm = np.empty((n_h + 1, n_span + 1, 3))

    for j in range(n_span + 1):
        # Global span index -> (panel, local w̄).  The last column of a
        # panel is the first column of the next, so clamp at the end.
        panel = min(j // n_w, params.N - 2)
        w = (j - panel * n_w) / n_w
        q1, q2 = q[panel], q[panel + 1]
        z0 = SPAN_SIGN * panel * params.D

        for i, h in enumerate(h_vals):
            pts[i, j] = surface_point(h, w, q1, q2, params.D, z0)
            nrm[i, j] = surface_normal_and_delta(h, w, q1, q2, params.D)[0]

    return pts, nrm


def build_fin_solid(
    params: FinParams,
    t: float = 0.0,
    n_h: int = 60,
    n_w: int = 24,
    thickness: float = 2.0e-3,
) -> tuple[np.ndarray, np.ndarray]:
    """Inflate the membrane into a closed triangulated solid.

    Parameters
    ----------
    params
        Fin parameters.
    t
        Time at which the ray angles are evaluated [s].  The STL is the
        *reference* configuration; motion is applied later by the
        ``pointDisplacement`` boundary condition.
    n_h, n_w
        Subdivisions along the height, and along each panel's span.
    thickness
        Membrane thickness [m].  Small compared with ``hmax``; it exists
        only so the surface can be closed.

    Returns
    -------
    vertices : ndarray, shape (n_vertices, 3)
    faces : ndarray, shape (n_faces, 3)
    """
    if thickness <= 0:
        raise ValueError("thickness must be positive")

    pts, nrm = _mid_surface(params, t, n_h, n_w)
    half = 0.5 * thickness

    top = pts + half * nrm      # +normal shell
    bot = pts - half * nrm      # −normal shell

    n_i, n_j = pts.shape[0], pts.shape[1]
    n_grid = n_i * n_j

    verts = np.concatenate([top.reshape(-1, 3), bot.reshape(-1, 3)], axis=0)

    def it(i: int, j: int) -> int:      # index into the top shell
        return i * n_j + j

    def ib(i: int, j: int) -> int:      # index into the bottom shell
        return n_grid + i * n_j + j

    faces: list[tuple[int, int, int]] = []

    # Two shells.  The bottom shell is wound the opposite way so its
    # facet normals point along −n, i.e. out of the solid.
    for i in range(n_i - 1):
        for j in range(n_j - 1):
            faces.append((it(i, j), it(i, j + 1), it(i + 1, j + 1)))
            faces.append((it(i, j), it(i + 1, j + 1), it(i + 1, j)))
            faces.append((ib(i, j), ib(i + 1, j + 1), ib(i, j + 1)))
            faces.append((ib(i, j), ib(i + 1, j), ib(i + 1, j + 1)))

    # Four side walls closing the rim.  Winding is chosen per edge so
    # that every wall's normal points away from the sheet.
    for j in range(n_j - 1):                       # h = hmin rim
        faces.append((it(0, j), it(0, j + 1), ib(0, j + 1)))
        faces.append((it(0, j), ib(0, j + 1), ib(0, j)))
    for j in range(n_j - 1):                       # h = hmax rim
        faces.append((it(n_i - 1, j), ib(n_i - 1, j), ib(n_i - 1, j + 1)))
        faces.append((it(n_i - 1, j), ib(n_i - 1, j + 1), it(n_i - 1, j + 1)))
    for i in range(n_i - 1):                       # w̄ = 0 rim (first ray)
        faces.append((it(i, 0), ib(i, 0), ib(i + 1, 0)))
        faces.append((it(i, 0), ib(i + 1, 0), it(i + 1, 0)))
    for i in range(n_i - 1):                       # w̄ = 1 rim (last ray)
        faces.append((it(i, n_j - 1), it(i + 1, n_j - 1), ib(i + 1, n_j - 1)))
        faces.append((it(i, n_j - 1), ib(i + 1, n_j - 1), ib(i, n_j - 1)))

    face_arr = np.asarray(faces, dtype=int)
    if not is_watertight(face_arr):
        raise RuntimeError(
            "fin surface is not watertight — every edge must be shared by "
            "exactly two triangles before it can be oriented"
        )
    return verts, orient_outward(verts, face_arr)


# ── Validity ────────────────────────────────────────────────────────

def is_watertight(faces: np.ndarray) -> bool:
    """True if every edge is shared by exactly two triangles.

    Necessary for ``snappyHexMesh`` to tell inside from outside.
    """
    from collections import Counter

    counts: Counter = Counter()
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            counts[(u, v) if u < v else (v, u)] += 1
    return all(n == 2 for n in counts.values())


# ── Writers ─────────────────────────────────────────────────────────

def _facet_normals(verts: np.ndarray, faces: np.ndarray) -> np.ndarray:
    a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    n = np.cross(b - a, c - a)
    mag = np.linalg.norm(n, axis=1, keepdims=True)
    return np.divide(n, mag, out=np.zeros_like(n), where=mag > 0)


def write_stl_ascii(
    path: str | Path, verts: np.ndarray, faces: np.ndarray, name: str = "fin",
) -> Path:
    """Write an ASCII STL.  Readable, but large — prefer binary."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    nrm = _facet_normals(verts, faces)

    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write(f"solid {name}\n")
        for f, n in zip(faces, nrm):
            fh.write(f"  facet normal {n[0]:.9e} {n[1]:.9e} {n[2]:.9e}\n")
            fh.write("    outer loop\n")
            for vi in f:
                v = verts[vi]
                fh.write(f"      vertex {v[0]:.9e} {v[1]:.9e} {v[2]:.9e}\n")
            fh.write("    endloop\n  endfacet\n")
        fh.write(f"endsolid {name}\n")
    return path


def write_stl_binary(
    path: str | Path, verts: np.ndarray, faces: np.ndarray, name: str = "fin",
) -> Path:
    """Write a binary STL."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    nrm = _facet_normals(verts, faces)

    with path.open("wb") as fh:
        fh.write(name.encode("ascii", "replace")[:80].ljust(80, b"\0"))
        fh.write(struct.pack("<I", len(faces)))
        for f, n in zip(faces, nrm):
            fh.write(struct.pack("<3f", *n))
            for vi in f:
                fh.write(struct.pack("<3f", *verts[vi]))
            fh.write(struct.pack("<H", 0))
    return path
