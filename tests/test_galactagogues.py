"""Tests for galactagogues and lactation induction.

These assert the pharmacological logic is sound: PK accumulation, prolactin
dose-response monotonicity and boundedness, and the protocol timeline shape.
Clinical accuracy is NOT asserted — the parameters are extrapolated from
postpartum cis literature and carry large uncertainty (documented in
docs/references.md).
"""

import numpy as np
import pytest

from hormone_sim import (
    GALACTAGOGUES,
    get_galactagogue,
    simulate_dosing_pk,
    prolactin_response,
    prolactin_curve,
    simulate_lactation,
    LactationProtocol,
    LactationModel,
    PROLACTIN_BASELINE_NG_ML,
)
from hormone_sim.galactagogues import list_galactagogues


# --- registry -----------------------------------------------------------------

def test_known_galactagogues_present():
    for name in ("metoclopramide", "domperidone", "sulpiride", "amisulpride"):
        g = get_galactagogue(name)
        assert g.half_life_h > 0
        assert 0 < g.bioavailability <= 1
        assert g.vd_l_kg > 0


def test_unknown_galactagogue_raises():
    with pytest.raises(KeyError):
        get_galactagogue("fenugreek")  # herbal: not a D2 antagonist, not in registry


def test_domperidone_not_us_approved():
    assert get_galactagogue("domperidone").approved_us is False
    assert get_galactagogue("metoclopramide").approved_us is True


def test_every_galactagogue_has_safety_note():
    for g in GALACTAGOGUES.values():
        assert g.safety_note.strip(), f"{g.name} missing safety note"


# --- PK -----------------------------------------------------------------------

def test_pk_accumulates_to_steady_state():
    g = get_galactagogue("domperidone")
    c = simulate_dosing_pk(g, 13.33, dosing_interval_h=8, duration_h=24 * 21, dt_h=0.1)
    # Steady state: mean of the last 7 days should exceed the single-dose peak.
    ss = c[-7 * 240:]
    assert ss.mean() > g.peak_conc_ng_ml(13.33), "accumulation should exceed single-dose peak"


def test_pk_zero_dose_is_zero():
    g = get_galactagogue("metoclopramide")
    c = simulate_dosing_pk(g, 0.0, dosing_interval_h=8, duration_h=24.0, dt_h=0.1)
    assert np.allclose(c, 0.0)


def test_pk_clearance_consistent_with_half_life():
    g = get_galactagogue("metoclopramide")
    cl = g.clearance_l_per_h(70.0)
    expected = np.log(2) * 2.5 * 70.0 / 5.5
    assert cl == pytest.approx(expected, rel=1e-6)


# --- prolactin pharmacodynamics ------------------------------------------------

def test_prolactin_baseline_at_zero_dose():
    g = get_galactagogue("domperidone")
    assert prolactin_response(g, 0.0) == pytest.approx(PROLACTIN_BASELINE_NG_ML)


def test_prolactin_monotonic_and_bounded():
    g = get_galactagogue("domperidone")
    doses = [0, 10, 30, 60, 90, 300, 3000]
    vals = [prolactin_response(g, d) for d in doses]
    assert all(vals[i] <= vals[i + 1] + 1e-9 for i in range(len(vals) - 1)), \
        "prolactin should be non-decreasing with dose"
    ceiling = PROLACTIN_BASELINE_NG_ML + g.prolactin_emax_ng_ml
    assert all(v <= ceiling + 1e-9 for v in vals), "should not exceed baseline+Emax"


def test_prolactin_curve_tracks_dose_response():
    # Time-varying curve's steady-state average should approximate the
    # dose-based steady-state value (rapid-equilibrium approximation).
    g = get_galactagogue("metoclopramide")
    dose = 30.0
    curve = prolactin_curve(g, dose, dosing_interval_h=8, duration_h=24 * 14, dt_h=0.1)
    ss_avg = float(curve[-7 * 240:].mean())
    target = prolactin_response(g, dose)
    assert ss_avg == pytest.approx(target, rel=0.15), \
        "curve average should match the dose-based steady state"


def test_prolactin_curve_never_below_baseline():
    g = get_galactagogue("domperidone")
    curve = prolactin_curve(g, 40.0, dosing_interval_h=8, duration_h=24 * 14, dt_h=0.1)
    assert np.all(curve >= PROLACTIN_BASELINE_NG_ML - 1e-9)


# --- lactation protocol ---------------------------------------------------------

def test_lactation_no_milk_during_prep():
    p = LactationProtocol(galactagogue="domperidone", galactagogue_dose_mg_day=40.0)
    r = simulate_lactation(p)
    prep_mask = r.phase == 0
    assert np.all(r.milk_ml_day[prep_mask] == 0.0), "no milk before the drop"


def test_lactation_milk_appears_after_drop():
    p = LactationProtocol(galactagogue="domperidone", galactagogue_dose_mg_day=40.0)
    r = simulate_lactation(p)
    lact_mask = r.phase == 1
    assert np.all(r.milk_ml_day[lact_mask] > 0.0), "milk after the drop"


def test_lactation_maturation_ramp_takes_weeks():
    """Milk should not be at full supply immediately after the drop."""
    p = LactationProtocol(galactagogue="domperidone", galactagogue_dose_mg_day=40.0)
    r = simulate_lactation(p)
    ss = r.steady_state_milk()
    drop = p.prep_weeks
    i_drop = int(np.argmax(r.weeks >= drop))
    i_1wk = int(round((drop + 1.0) / 0.25))
    i_8wk = int(round((drop + 8.0) / 0.25))
    assert r.milk_ml_day[i_1wk] < r.milk_ml_day[i_8wk] <= r.milk_ml_day.max()
    # Early supply is a fraction of steady state, not the whole thing.
    assert r.milk_ml_day[i_drop] < 0.25 * ss


def test_lactation_more_pumping_more_milk():
    low = simulate_lactation(LactationProtocol(pumping_sessions_per_day=4.0))
    high = simulate_lactation(LactationProtocol(pumping_sessions_per_day=12.0))
    assert high.steady_state_milk() > low.steady_state_milk()


def test_lactation_lower_dose_lower_milk():
    low = simulate_lactation(LactationProtocol(galactagogue_dose_mg_day=10.0))
    high = simulate_lactation(LactationProtocol(galactagogue_dose_mg_day=60.0))
    assert high.steady_state_milk() >= low.steady_state_milk()


def test_lactation_priming_grows_during_prep():
    p = LactationProtocol(prep_weeks=30)
    r = simulate_lactation(p)
    prep = r.priming[r.phase == 0]
    assert np.all(np.diff(prep) >= 0), "priming should be non-decreasing during prep"
    assert prep[-1] > 0.5, "30 weeks of prep should prime substantially"


def test_lactation_adequate_prep_needed():
    """Minimal preparation should give a much lower ceiling than full prep."""
    short = simulate_lactation(LactationProtocol(prep_weeks=4))
    full = simulate_lactation(LactationProtocol(prep_weeks=30))
    assert full.steady_state_milk() > short.steady_state_milk()


def test_lactation_unknown_galactagogue_raises():
    with pytest.raises(KeyError):
        simulate_lactation(LactationProtocol(galactagogue="blessed_thistle"))


def test_list_galactagogues_runs():
    txt = list_galactagogues()
    assert "domperidone" in txt
    assert "Safety" in txt
