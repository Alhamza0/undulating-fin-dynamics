"""The comparison layer decides whether the CFD campaign succeeded, so
its windowing and metrics need to be right before any result rests on
them.
"""

import numpy as np
import pytest

from undulating_fin.params import FinParams
from openfoam.compare import (
    COMPONENTS,
    ForceSeries,
    analytical_forces,
    is_periodic,
    last_cycle,
    load_openfoam_forces,
    metrics,
)


P = FinParams(N=2, D=0.09, hmin=0.034, hmax=0.174,
              A=np.deg2rad(30), f=1.0, Phi=0.0)


def _write_force_dat(d, t, F, time_dir="0"):
    out = d / "postProcessing/forces" / time_dir
    out.mkdir(parents=True)
    rows = np.column_stack([t, F, F, np.zeros_like(F)])
    with (out / "force.dat").open("w", newline="\n") as fh:
        fh.write("# Force\n# CofR : (0 0 0)\n#\n")
        fh.write("# Time\ttotal_x total_y total_z\tpressure viscous\n")
        for r in rows:
            fh.write("\t".join("%.8e" % v for v in r) + "\n")
    return d


# ── Loading ─────────────────────────────────────────────────────────

def test_reads_force_dat(tmp_path):
    t = np.linspace(0, 1, 11)
    F = np.column_stack([t, 2 * t, 3 * t])
    s = load_openfoam_forces(_write_force_dat(tmp_path, t, F))
    assert s.t == pytest.approx(t)
    assert s.F == pytest.approx(F)


def test_concatenates_restarts_and_drops_duplicates(tmp_path):
    """A restarted run writes overlapping time directories; the later
    values must win rather than appearing twice."""
    t1 = np.linspace(0.0, 1.0, 11)
    t2 = np.linspace(1.0, 2.0, 11)          # 1.0 appears in both
    _write_force_dat(tmp_path, t1, np.column_stack([t1] * 3), "0")
    _write_force_dat(tmp_path, t2, np.column_stack([t2 * 10] * 3), "1")

    s = load_openfoam_forces(tmp_path)
    assert len(s) == 21
    assert np.all(np.diff(s.t) > 0)
    assert s.F[s.t == 1.0][0, 0] == pytest.approx(10.0)   # the restart's value


def test_missing_output_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_openfoam_forces(tmp_path)


def test_shape_mismatch_rejected():
    with pytest.raises(ValueError):
        ForceSeries(np.zeros(5), np.zeros((4, 3)))


# ── Windowing ───────────────────────────────────────────────────────

def test_last_cycle_takes_exactly_one_period():
    t = np.linspace(0, 3, 3001)
    s = ForceSeries(t, np.zeros((t.size, 3)))
    w = last_cycle(s, f=1.0, n_cycles=1)
    assert w.t[0] == pytest.approx(2.0)
    assert w.t[-1] == pytest.approx(3.0)


def test_last_cycle_refuses_a_run_that_is_too_short():
    t = np.linspace(0, 0.5, 51)
    s = ForceSeries(t, np.zeros((t.size, 3)))
    with pytest.raises(ValueError, match="too short"):
        last_cycle(s, f=1.0, n_cycles=1)


def test_transient_is_detected_as_not_periodic():
    """A decaying start-up transient must not pass as a periodic state.

    This is the guard that stops a comparison measuring how the fluid
    was set in motion rather than how the fin drives it.
    """
    t = np.linspace(0, 3, 3001)
    clean = np.sin(2 * np.pi * t)
    transient = clean + 5.0 * np.exp(-t)
    assert is_periodic(ForceSeries(t, np.column_stack([clean] * 3)), 1.0)
    assert not is_periodic(ForceSeries(t, np.column_stack([transient] * 3)), 1.0)


# ── Metrics ─────────────────────────────────────────────────────────

def test_identical_series_score_perfectly():
    t = np.linspace(0, 1, 101)
    F = np.column_stack([np.sin(2 * np.pi * t), np.cos(2 * np.pi * t), t])
    m = metrics(ForceSeries(t, F), ForceSeries(t, F))
    for c in COMPONENTS:
        assert m[c]["rmse_N"] == pytest.approx(0.0, abs=1e-15)
        assert m[c]["nrmse_pct"] == pytest.approx(0.0, abs=1e-12)
        assert m[c]["correlation"] == pytest.approx(1.0)


def test_sign_flip_is_reported_as_disagreement():
    """An inverted model must score correlation -1, not +1.

    Both sides report the force on the fin, so a sign difference is a
    real disagreement and must never be silently absorbed.
    """
    t = np.linspace(0, 1, 101)
    F = np.column_stack([np.sin(2 * np.pi * t)] * 3)
    m = metrics(ForceSeries(t, F), ForceSeries(t, -F))
    for c in COMPONENTS:
        assert m[c]["correlation"] == pytest.approx(-1.0)
        assert m[c]["rmse_N"] > 0


def test_nrmse_is_normalised_by_the_cfd_range():
    t = np.linspace(0, 1, 101)
    cfd = np.column_stack([np.sin(2 * np.pi * t)] * 3)      # range 2
    model = cfd + 0.2
    m = metrics(ForceSeries(t, cfd), ForceSeries(t, model))
    for c in COMPONENTS:
        assert m[c]["nrmse_pct"] == pytest.approx(10.0, rel=1e-9)


def test_constant_series_gives_nan_correlation_not_a_crash():
    t = np.linspace(0, 1, 51)
    z = np.zeros((t.size, 3))
    m = metrics(ForceSeries(t, z), ForceSeries(t, z))
    assert np.isnan(m["surge"]["correlation"])


def test_mismatched_time_bases_rejected():
    a = ForceSeries(np.linspace(0, 1, 10), np.zeros((10, 3)))
    b = ForceSeries(np.linspace(0, 2, 10), np.zeros((10, 3)))
    with pytest.raises(ValueError, match="time base"):
        metrics(a, b)


# ── Analytical sampling ─────────────────────────────────────────────

def test_analytical_matches_the_requested_time_base():
    t = np.linspace(2.0, 3.0, 41)
    s = analytical_forces(P, t, "morison")
    assert s.t == pytest.approx(t)
    assert np.all(np.isfinite(s.F))


def test_surge_vanishes_at_zero_phase_shift():
    """With Phi = 0 the rays move together, so there is no travelling
    wave and no net thrust - the physical check the probe already
    passed on the CFD side."""
    t = np.linspace(2.0, 3.0, 101)
    for method in ("rft", "morison"):
        s = analytical_forces(P, t, method)
        assert np.abs(s.F[:, 2]).max() < 1e-9, method


# ── The added-mass trap ─────────────────────────────────────────────

def test_morison_requires_accelerations():
    """Omitting ddq must fail loudly, not silently become RFT.

    compute_total_forces used to substitute zeros, which removes the
    added-mass term entirely - so 'morison' returned RFT's answer under
    a different name. The paper's headline result is that the two models
    differ, so a silent collapse is the worst possible failure.
    """
    from undulating_fin.hydrodynamics import compute_total_forces

    q = np.array([0.3, 0.1])
    dq = np.array([1.0, -0.5])
    for method in ("morison", "lighthill"):
        with pytest.raises(ValueError, match="needs joint accelerations"):
            compute_total_forces(P, q, dq, method=method)

    # RFT genuinely does not use them, so it must still work.
    assert compute_total_forces(P, q, dq, method="rft").shape == (3,)


def test_morison_and_rft_actually_differ():
    """With accelerations supplied, the two models must separate.

    Guards the fix above: if this ever collapses to equality again, the
    added-mass term has gone missing.
    """
    t = np.linspace(2.0, 3.0, 101)
    p = FinParams(N=2, D=0.09, hmin=0.034, hmax=0.174, A=np.deg2rad(30),
                  f=1.0, Phi=np.deg2rad(40), Cn=2.7, Cm=1.0, rho=1000.0)
    rft = analytical_forces(p, t, "rft").F
    mor = analytical_forces(p, t, "morison").F
    assert np.abs(rft - mor).max() > 1e-3
