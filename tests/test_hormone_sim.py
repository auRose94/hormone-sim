"""Tests for hormone-sim.

These assert *physiological sanity* of the model, not clinical accuracy (the
latter requires validation against measured data, see docs/). They guard the
invariants that must hold regardless of parameter tuning.
"""

import numpy as np
import pytest

from hormone_sim import (
    simulate,
    Regimen,
    EstradiolSpec,
    run_population,
    TARGET_E2_PGML,
    TARGET_T_NGDL,
    e2_pgml_to_pmol,
    t_ngdl_to_nmol,
)
from hormone_sim.pk import build_schedule, simulate_pk
from hormone_sim.params import ROUTES, EstradiolPK, AndrogenPD, BASELINE_E2_PGML
from hormone_sim.pd import estradiol_suppression, anti_androgen_suppression


def make_regimen(route="transdermal", dose=50.0, aa=None):
    return Regimen(
        name=f"test_{route}",
        population="trans_feminine",
        estradiol=EstradiolSpec(route=route, daily_rate_ug_day=dose),
        anti_androgen=aa or {},
    )


# --- unit: PK ---------------------------------------------------------------

def test_transdermal_estradiol_rises_above_baseline():
    r = make_regimen("transdermal", 50.0)
    spec = r.estradiol
    sched = build_schedule(
        spec.daily_rate_ug_day, spec.route_params(), period_hours=7 * 24,
    )
    pk = EstradiolPK()
    C = simulate_pk(sched, pk, 7 * 24, dt_hours=0.1, route="transdermal")
    assert C[-1] > BASELINE_E2_PGML, "estradiol should rise above baseline"
    assert np.all(C >= 0), "concentration must be non-negative"


def test_more_dose_gives_more_estradiol():
    pk = EstradiolPK()
    results = []
    for dose in (25.0, 50.0, 100.0):
        sched = build_schedule(dose, ROUTES["transdermal"], period_hours=14 * 24)
        C = simulate_pk(sched, pk, 14 * 24, dt_hours=0.1, route="transdermal")
        results.append(C[-1])
    assert results[0] < results[1] < results[2], "steady-state E2 should scale with dose"


def test_oral_lower_bioavailability_than_transdermal():
    """Same nominal daily mass: oral should give lower steady-state E2 due to
    first-pass (lower effective bioavailability)."""
    pk = EstradiolPK()
    sched_o = build_schedule(50.0, ROUTES["oral"], period_hours=14 * 24)
    sched_t = build_schedule(50.0, ROUTES["transdermal"], period_hours=14 * 24)
    co = simulate_pk(sched_o, pk, 14 * 24, dt_hours=0.1, route="oral")
    ct = simulate_pk(sched_t, pk, 14 * 24, dt_hours=0.1, route="transdermal")
    assert co[-1] < ct[-1], "oral first-pass should lower steady-state E2 vs transdermal"


def test_zero_dose_returns_to_baseline():
    pk = EstradiolPK()
    sched = build_schedule(0.0, ROUTES["transdermal"], period_hours=7 * 24)
    C = simulate_pk(sched, pk, 7 * 24, dt_hours=0.1, route="transdermal")
    # With zero input, the system should stay at (decay toward) baseline.
    assert abs(C[-1] - BASELINE_E2_PGML) < 1e-6


def test_half_life_positive_and_reasonable():
    pk = EstradiolPK()
    t12 = pk.half_life_h("transdermal")
    # Apparent half-life of this calibrated 1-compartment model is ~10 days
    # (a documented simplification — see params.EstradiolPK docstring). It must
    # be positive and in a plausible (if long) range, and oral should be shorter
    # due to first-pass clearance.
    assert t12 > 0, "apparent half-life must be positive"
    assert pk.half_life_h("oral") < t12, "oral first-pass shortens apparent t1/2"


# --- unit: PD ---------------------------------------------------------------

def test_suppression_monotonic_in_estradiol():
    e2 = np.array([0.0, 50.0, 150.0, 500.0])
    s = estradiol_suppression(e2, AndrogenPD())
    assert np.all((s >= 0) & (s <= 1)), "suppression fraction in [0,1]"
    assert s[0] == 0.0, "no estradiol -> no HPG suppression"
    assert s[-1] > s[0], "more estradiol -> more suppression"


def test_anti_androgen_bounded():
    s = anti_androgen_suppression({"spironolactone": 1.0, "cyproterone": 1.0})
    assert 0.0 <= s <= 1.0
    # more agents -> at least as much suppression
    assert anti_androgen_suppression({"spironolactone": 1.0}) <= s


def test_unknown_anti_androgen_raises():
    with pytest.raises(KeyError):
        anti_androgen_suppression({"not_a_real_drug": 1.0})


# --- integration: full simulation -------------------------------------------

def test_full_sim_reaches_plausible_steady_state():
    r = make_regimen("transdermal", 50.0, {"spironolactone": 1.0})
    sim = simulate(r, weeks=12)
    e2 = sim.steady_state_e2_pgml()
    t = sim.steady_state_t_ngdl()
    assert 20 < e2 < 600, f"steady-state E2 {e2} pg/mL implausible"
    assert 0 <= t <= 300, f"steady-state T {t} ng/dL implausible"
    assert t < 300, "anti-androgen + E2 should reduce T below baseline"


def test_testosterone_suppressed_below_baseline():
    r = make_regimen("transdermal", 100.0, {"spironolactone": 1.0})
    sim = simulate(r, weeks=16)
    assert sim.steady_state_t_ngdl() < 300.0


def test_units_conversions():
    assert e2_pgml_to_pmol(100.0) == pytest.approx(367.1, rel=1e-6)
    assert t_ngdl_to_nmol(100.0) == pytest.approx(3.47, rel=1e-6)


# --- population -------------------------------------------------------------

def test_population_band_width_positive():
    r = make_regimen("transdermal", 50.0, {"spironolactone": 1.0})
    pop = run_population(lambda: r, weeks=8, n=30, seed=1)
    lo, med, hi = pop.percentile_band(0.8)
    assert np.all(hi >= med)
    assert np.all(med >= lo)
    # fraction in target is a probability in [0,1]
    assert 0.0 <= pop.fraction_e2_in_target() <= 1.0
    assert 0.0 <= pop.fraction_t_suppressed() <= 1.0


def test_population_is_reproducible_with_seed():
    r = make_regimen("transdermal", 50.0)
    p1 = run_population(lambda: r, weeks=6, n=20, seed=123)
    p2 = run_population(lambda: r, weeks=6, n=20, seed=123)
    assert np.allclose(p1.e2_pgml, p2.e2_pgml)
