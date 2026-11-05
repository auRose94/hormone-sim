"""Tests for equieffective comparison, progesterone phase, sensitivity bands."""

import numpy as np
import pytest

from hormone_sim import (
    get_galactagogue,
    equieffective_dose,
    compare_galactagogues,
    comparison_table,
    tradeoff_notes,
    LactationProtocol,
    LactationModel,
    simulate_lactation,
    SensitivitySpec,
    simulate_lactation_population,
    PROLACTIN_BASELINE_NG_ML,
)


# --- equieffective dose / comparison ------------------------------------------

def test_equieffective_dose_reaches_target():
    for name in ("domperidone", "metoclopramide", "sulpiride"):
        g = get_galactagogue(name)
        d = equieffective_dose(g, 130.0)
        assert d > 0
        assert PROLACTIN_BASELINE_NG_ML < 130.0  # sanity
        # Reconstruct: prolactin at that dose should equal the target.
        rise = 130.0 - PROLACTIN_BASELINE_NG_ML
        back = g.prolactin_emax_ng_ml * d / (g.prolactin_ed50_mg_day + d)
        assert back == pytest.approx(rise, rel=1e-6)


def test_equieffective_dose_above_ceiling_raises():
    g = get_galactagogue("metoclopramide")  # lowest Emax (150)
    ceiling = PROLACTIN_BASELINE_NG_ML + g.prolactin_emax_ng_ml
    with pytest.raises(ValueError):
        equieffective_dose(g, ceiling + 50.0)


def test_compare_same_supply_at_equieffective_prolactin():
    """The model's core claim: at equal prolactin, equal supply — differences
    are PK/safety, not milk output."""
    rows = compare_galactagogues(["domperidone", "metoclopramide"],
                                 target_prolactin_ng_ml=130.0)
    assert len(rows) == 2
    assert rows[0].ss_milk == pytest.approx(rows[1].ss_milk, rel=0.02)
    # Doses differ (different potency).
    assert rows[0].dose_mg_day != rows[1].dose_mg_day


def test_compare_requires_more_dose_for_weaker_agent():
    d_domp = equieffective_dose(get_galactagogue("domperidone"), 130.0)
    d_meto = equieffective_dose(get_galactagogue("metoclopramide"), 130.0)
    # metoclopramide has lower Emax and higher ED50 -> needs more dose here.
    assert d_meto > d_domp


def test_comparison_table_and_notes_run():
    rows = compare_galactagogues(["domperidone", "metoclopramide"])
    t = comparison_table(rows)
    assert "domperidone" in t and "metoclopramide" in t
    n = tradeoff_notes(rows)
    assert "steady state" in n


# --- progesterone phase (two-compartment priming) -------------------------------

def test_la_requires_progesterone():
    """Without a P phase, lobuloalveolar development never starts (model
    assumption: LA is P-driven)."""
    r = simulate_lactation(LactationProtocol(
        prep_weeks=30, progesterone_start_weeks=30))  # P never active
    assert np.all(r.lobuloalveolar == 0.0)
    assert np.all(r.milk_ml_day == 0.0)
    assert r.ductal.max() > 0.5  # ducts still grew on E2


def test_longer_p_phase_more_capacity():
    caps = {}
    for p_start in (26.0, 18.0, 10.0):
        r = simulate_lactation(LactationProtocol(
            prep_weeks=30, progesterone_start_weeks=p_start))
        caps[p_start] = r.priming[int(round(30 / 0.25))]
    assert caps[10.0] > caps[18.0] > caps[26.0]


def test_p_window_shown_in_result():
    r = simulate_lactation(LactationProtocol(
        prep_weeks=30, progesterone_start_weeks=18))
    i_before = int(round(10 / 0.25))
    i_during = int(round(25 / 0.25))
    assert r.lobuloalveolar[i_before] == 0.0
    assert r.lobuloalveolar[i_during] > 0.0


def test_ductal_grows_with_e2_level():
    hi = simulate_lactation(LactationProtocol(prep_weeks=20, prep_e2_pgml=300))
    lo = simulate_lactation(LactationProtocol(prep_weeks=20, prep_e2_pgml=75))
    i = int(round(20 / 0.25))
    assert hi.ductal[i] > lo.ductal[i]


def test_capacity_geometric_mean():
    r = simulate_lactation(LactationProtocol(
        prep_weeks=30, progesterone_start_weeks=18))
    i = int(round(30 / 0.25))
    expected = np.sqrt(r.ductal[i] * r.lobuloalveolar[i])
    assert r.priming[i] == pytest.approx(expected, rel=1e-6)


def test_default_p_start_is_prep_minus_12():
    p = LactationProtocol(prep_weeks=30)
    assert p.p_start() == 18.0
    p2 = LactationProtocol(prep_weeks=8)  # shorter than the 12-wk default
    assert p2.p_start() == 0.0


# --- sensitivity bands ----------------------------------------------------------

def test_population_bands_reproducible_with_seed():
    p = LactationProtocol(prep_weeks=20)
    a = simulate_lactation_population(p, n=20, seed=7)
    b = simulate_lactation_population(p, n=20, seed=7)
    assert np.allclose(a.milk, b.milk)


def test_population_bands_contain_median_and_spread():
    p = LactationProtocol(prep_weeks=30)
    pop = simulate_lactation_population(p, n=100, seed=42)
    lo, med, hi = pop.band(0.8)
    assert np.all(lo <= med) and np.all(med <= hi)
    # Spread should be real (not degenerate).
    assert (hi - lo)[-1] > 10.0


def test_population_steady_state_ordering():
    p = LactationProtocol(prep_weeks=30)
    pop = simulate_lactation_population(p, n=100, seed=42)
    p10, p50, p90 = pop.steady_state()
    assert p10 < p50 < p90


def test_population_default_median_near_deterministic():
    """With modest CVs, the population median should sit near the
    deterministic run (lognormal sampling has median-preserving center)."""
    p = LactationProtocol(prep_weeks=30, galactagogue="domperidone",
                          galactagogue_dose_mg_day=40.0)
    det = simulate_lactation(p)
    pop = simulate_lactation_population(p, n=100, seed=42)
    _, med, _ = pop.steady_state()
    # Geometric-mean-of-CVs is modest; allow generous 25% tolerance because
    # milk is a product of several lognormal terms (compound variance).
    assert med == pytest.approx(det.steady_state_milk(), rel=0.25)


def test_maturation_tau_has_widest_band():
    spec = SensitivitySpec()
    assert spec.cvs["maturation_tau_weeks"] > max(
        v for k, v in spec.cvs.items() if k != "maturation_tau_weeks")


def test_population_report_runs():
    pop = simulate_lactation_population(LactationProtocol(prep_weeks=20), n=30)
    txt = pop.report()
    assert "sensitivity bands" in txt
    assert "no controlled trials" in txt
