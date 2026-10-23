#!/usr/bin/env python3
"""
Compute fin boundary displacement at each time step.

Usage:
    python fin_motion.py <time> <points_file> <output_file>

The script reads the reference fin boundary points, applies the
traveling-wave kinematics from FinParams, and writes the
displacement vectors.
"""

import sys
import numpy as np

# Fin parameters (frozen from case generation)
N      = 3
A      = 0.5235987755982988         # amplitude [rad]
f      = 2.0         # frequency [Hz]
Phi    = 0.8726646259971648       # inter-ray phase shift [rad]
X      = 0.0         # angular offset [rad]
D      = 0.1         # inter-ray distance [m]
hmax   = 0.175      # max membrane height [m]
hmin   = 0.0      # min membrane height [m]


def ray_angle(i: int, t: float) -> float:
    """Joint angle for ray *i* (0-indexed) at time *t*."""
    return A * np.sin(2 * np.pi * f * t - i * Phi) + X


def membrane_displacement(x_ref: float, y_ref: float, t: float):
    """Compute (dx, dy) displacement for a mesh point at (x_ref, y_ref).

    Determines which panel the point belongs to from its x-coordinate,
    computes the two bounding ray angles, and interpolates the
    displacement using the normalized span position w_bar.
    """
    # Identify panel index
    panel = int(x_ref / D)
    panel = max(0, min(panel, N - 2))
    q1 = ray_angle(panel, t)
    q2 = ray_angle(panel + 1, t)

    # Normalized span
    x_local = x_ref - panel * D
    w_bar = x_local / D if D > 0 else 0.0
    w_bar = max(0.0, min(1.0, w_bar))

    # Height (distance from joint axis)
    h = abs(y_ref)
    h = max(hmin, min(h, hmax))
    sign_y = 1.0 if y_ref >= 0 else -1.0

    # Displacement from rotating rays
    # Ray at panel rotates by q1, ray at panel+1 by q2
    # Interpolated angle at this span position
    q_interp = q1 * (1.0 - w_bar) + q2 * w_bar

    # Displacement of this point (rotation about joint)
    # Reference position: y = h, x = x_local relative to joint
    dx = 0.0  # no axial displacement (rays rotate in y)
    dy = sign_y * h * (np.cos(q_interp) - 1.0)  # transverse displacement

    # Add the tangential component
    dx += -sign_y * h * np.sin(q_interp)

    return dx, dy


def main():
    if len(sys.argv) < 4:
        print("Usage: fin_motion.py <time> <points_file> <output_file>")
        sys.exit(1)

    t = float(sys.argv[1])
    points_file = sys.argv[2]
    output_file = sys.argv[3]

    # Read reference points (OpenFOAM format)
    points = []
    with open(points_file, "r") as fh:
        in_block = False
        for line in fh:
            line = line.strip()
            if line == "(":
                in_block = True
                continue
            if line == ")":
                break
            if in_block and line.startswith("("):
                coords = line.strip("()").split()
                points.append([float(c) for c in coords])
    points = np.array(points)

    # Compute displacements
    displacements = np.zeros_like(points)
    for idx in range(len(points)):
        dx, dy = membrane_displacement(points[idx, 0], points[idx, 1], t)
        displacements[idx, 0] = dx
        displacements[idx, 1] = dy
        # z = 0 (2-D)

    # Write OpenFOAM-format vector field
    with open(output_file, "w") as fh:
        fh.write(str(len(displacements)) + "\n")
        fh.write("(\n")
        for d in displacements:
            fh.write(f"({d[0]:.10e} {d[1]:.10e} {d[2]:.10e})\n")
        fh.write(")\n")

    print(f"Wrote {len(displacements)} displacement vectors for t = {t}")


if __name__ == "__main__":
    main()
