"""
Parse OpenFOAM postProcessing output into numpy arrays.

Handles the standard ``forces`` function object output
(``postProcessing/forces/0/force.dat``) as well as custom
coefficient files.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import NamedTuple

import numpy as np


class ForceData(NamedTuple):
    """Parsed force time-history."""
    time: np.ndarray          # (nT,)
    pressure: np.ndarray      # (nT, 3) — pressure contribution [Fx, Fy, Fz]
    viscous: np.ndarray       # (nT, 3) — viscous contribution [Fx, Fy, Fz]
    total: np.ndarray         # (nT, 3) — total = pressure + viscous


# Regex for OpenFOAM force.dat lines:
#   <time> (( px py pz ) ( vx vy vz ) ( px py pz ))
_FORCE_LINE = re.compile(
    r"^\s*"
    r"([\d.eE+\-]+)\s+"               # time
    r"\(\(\s*"
    r"([\d.eE+\-]+)\s+"               # pressure x
    r"([\d.eE+\-]+)\s+"               # pressure y
    r"([\d.eE+\-]+)\s*\)\s*\(\s*"     # pressure z
    r"([\d.eE+\-]+)\s+"               # viscous x
    r"([\d.eE+\-]+)\s+"               # viscous y
    r"([\d.eE+\-]+)\s*\)\s*\(\s*"     # viscous z
    r"([\d.eE+\-]+)\s+"               # porous x
    r"([\d.eE+\-]+)\s+"               # porous y
    r"([\d.eE+\-]+)\s*\)\)"           # porous z
)

_FORCE_LINE_V2 = re.compile(
    r"^\s*"
    r"([\d.eE+\-]+)\s+"               # time
    r"\(\s*"
    r"\(\s*([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\)\s*"  # pressure
    r"\(\s*([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\)\s*"  # viscous
    r"\(\s*([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\)\s*"  # porous
    r"\)"
)


def parse_forces(force_file: str | Path) -> ForceData:
    """Parse an OpenFOAM ``force.dat`` / ``forces.dat`` file.

    Parameters
    ----------
    force_file : path-like
        Path to the force data file, typically at
        ``<case>/postProcessing/forces/0/force.dat``.

    Returns
    -------
    ForceData
        Named tuple with time, pressure, viscous, and total arrays.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If no valid data lines are found.
    """
    path = Path(force_file)
    if not path.is_file():
        raise FileNotFoundError(f"Force file not found: {path}")

    times, pres, visc = [], [], []

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            m = _FORCE_LINE.match(line) or _FORCE_LINE_V2.match(line)
            if m is None:
                continue

            vals = [float(m.group(i)) for i in range(1, 11)]
            times.append(vals[0])
            pres.append(vals[1:4])
            visc.append(vals[4:7])
            # porous contribution (vals[7:10]) ignored for now

    if not times:
        raise ValueError(f"No valid force data found in {path}")

    t = np.array(times)
    p = np.array(pres)
    v = np.array(visc)
    return ForceData(time=t, pressure=p, viscous=v, total=p + v)


def find_force_file(case_dir: str | Path) -> Path:
    """Locate the force.dat file inside a case's postProcessing directory.

    Searches common paths:
    - postProcessing/forces/0/force.dat
    - postProcessing/forces/0/forces.dat
    - postProcessing/forceCoeffs/0/coefficient.dat
    """
    case = Path(case_dir)
    candidates = [
        case / "postProcessing" / "forces" / "0" / "force.dat",
        case / "postProcessing" / "forces" / "0" / "forces.dat",
        case / "postProcessing" / "forceCoeffs" / "0" / "coefficient.dat",
    ]

    # Also check for time-stamped subdirectories
    pp_forces = case / "postProcessing" / "forces"
    if pp_forces.is_dir():
        for subdir in sorted(pp_forces.iterdir()):
            if subdir.is_dir():
                for name in ("force.dat", "forces.dat"):
                    candidates.append(subdir / name)

    for c in candidates:
        if c.is_file():
            return c

    raise FileNotFoundError(
        f"No force data found in {case / 'postProcessing'}. "
        "Ensure the simulation has completed and force function objects are configured."
    )


def load_forces(case_dir: str | Path) -> ForceData:
    """Convenience: find and parse forces from a case directory."""
    return parse_forces(find_force_file(case_dir))


def resample_to_analytical(
    of_data: ForceData,
    time_analytical: np.ndarray,
) -> np.ndarray:
    """Linearly interpolate OpenFOAM forces onto the analytical time grid.

    Parameters
    ----------
    of_data : ForceData
        Parsed OpenFOAM results.
    time_analytical : 1-D array
        Target time vector from the analytical model.

    Returns
    -------
    F_of_resampled : ndarray, shape (len(time_analytical), 3)
        Interpolated [Fx, Fy, Fz] forces.
    """
    F = np.zeros((len(time_analytical), 3))
    for k in range(3):
        F[:, k] = np.interp(time_analytical, of_data.time, of_data.total[:, k])
    return F


def list_completed_cases(base_dir: str | Path) -> list[Path]:
    """Find all case directories under *base_dir* that have force results."""
    base = Path(base_dir)
    cases = []
    if not base.is_dir():
        return cases
    for d in sorted(base.iterdir()):
        if not d.is_dir():
            continue
        try:
            find_force_file(d)
            cases.append(d)
        except FileNotFoundError:
            continue
    return cases
