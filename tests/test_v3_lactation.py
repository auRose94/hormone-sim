"""Tests for progestogen PK, taper/brake behaviors, end-to-end bands, and
participant personalization."""

import numpy as np
import pytest

from hormone_sim import (
    get_progestogen,
    custom_progestogen,
    LactationProtocol,
    LactationModel,
    simulate_lactation,
    simulate_lactation_population,
    SensitivitySpec,
    personalize,
    posterior_predictive,
    Observation,
)


# --- progestogen PK -----------------------------------------------------------

def test_micronized_p_wild_peaks():
    """The headline PK fact: q24h dosing with t1/2 ~6 h gives ~3x peak/mean."""
    p = get_progestogen("micronized_progesterone")
    ratio = p.peak_to_mean_ratio(200.0)
    assert 2.0 < ratio < 5.0


def test_micronized_p_mean_level_matches_literature_scale():
    """200 mg/day oral micronized P -> mean ~10 ng/mL (published range ~5-20)."""
    p = get_progestogen("micronized_progesterone")
    assert 5.0 < p.mean_conc_ng_ml(200.0) < 20.0


def test_p_levels_scale_linearly_with_dose():
    p = get_progestogen("micronized_progesterone")
    assert p.mean_conc_ng_ml(300.0) == pytest.approx(
        3 * p.mean_conc_ng_ml(100.0), rel=1e-6)


def test_unknown_progestogen_raises():
    with pytest.raises(KeyError):
        get_progestogen("medroxyprogesterone")  # not characterized in registry


def test_custom_progestogen():
    c = custom_progestogen("test_synthetic", half_life_h=24.0,
                           bioavailability=1.0, vd_l_kg=4.0,
                           typical_dose_mg_day=10.0)
    assert c.mean_conc_ng_ml(10.0) > 0
    # Long t1/2 -> far less peaky than micronized P.
    assert c.peak_to_mean_ratio(10.0) < get_progestogen(
        "micronized_progesterone").peak_to_mean_ratio(200.0)


# --- P dose-response on lobuloalveolar development -----------------------------

def test_p_dose_drives_la_development():
    low = simulate_lactation(LactationProtocol(progesterone_dose_mg_day=100.0))
    high = simulate_lactation(LactationProtocol(progesterone_dose_mg_day=300.0))
    i = int(round(30 / 0.25))
    assert high.lobuloalveolar[i] > low.lobuloalveolar[i]
    assert high.steady_state_milk() >= low.steady_state_milk()


def test_p_free_protocol_still_zero_milk():
    r = simulate_lactation(LactationProtocol(
        prep_weeks=30, progesterone_start_weeks=30))
    assert np.all(r.milk_ml_day == 0.0)


def test_p_withdrawal_gate_unit():
    m = LactationModel()
    # Prep mean 50 ng/mL: withdrawal requires falling below 20% (10 ng/mL).
    assert not m.p_withdrawn(20.0, 50.0)
    assert m.p_withdrawn(5.0, 50.0)


# --- taper vs hard drop ----------------------------------------------------------

def test_taper_delays_onset_but_not_final_supply():
    hard = simulate_lactation(LactationProtocol(taper_weeks=0.0))
    soft = simulate_lactation(LactationProtocol(taper_weeks=4.0))
    wk_hard = hard.weeks_to_first_milk()
    wk_soft = soft.weeks_to_first_milk()
    assert wk_soft is not None and wk_hard is not None
    assert wk_soft > wk_hard, "tapering P delays lactogenesis II"
    # Final supply unchanged once fully matured and P withdrawn either way.
    assert soft.steady_state_milk() == pytest.approx(
        hard.steady_state_milk(), rel=0.1)


def test_taper_trajectories_are_linear_down():
    r = simulate_lactation(LactationProtocol(taper_weeks=4.0))
    # At drop+2 wk (mid-taper), P should be ~half its prep mean.
    prep_mean = r.p_level_ng_ml[int(round(25 / 0.25))]
    mid = r.p_level_ng_ml[int(round(32 / 0.25))]  # drop + 2 wk
    assert 0.3 * prep_mean < mid < 0.7 * prep_mean


# --- maintained E2 brake -----------------------------------------------------------

def test_post_drop_e2_brake_reduces_supply():
    none = simulate_lactation(LactationProtocol(post_drop_e2_pgml=0.0))
    kept = simulate_lactation(LactationProtocol(post_drop_e2_pgml=150.0))
    assert kept.steady_state_milk() < none.steady_state_milk()
    # And the brake is meaningful but not annihilating at HRT-level E2.
    ratio = kept.steady_state_milk() / none.steady_state_milk()
    assert 0.4 < ratio < 0.95


def test_no_brake_at_zero_e2():
    m = LactationModel()
    assert m.e2_brake(0.0) == pytest.approx(1.0)


# --- end-to-end bands --------------------------------------------------------------

def test_end_to_end_bands_wider_than_gland_only():
    p = LactationProtocol(galactagogue="metoclopramide",
                          galactagogue_dose_mg_day=15)
    both = simulate_lactation_population(p, n=120, seed=42,
                                         include_drug_uncertainty=True)
    gland = simulate_lactation_population(p, n=120, seed=42,
                                          include_drug_uncertainty=False)
    w_both = both.steady_state()[2] - both.steady_state()[0]
    w_gland = gland.steady_state()[2] - gland.steady_state()[0]
    assert w_both > w_gland, "drug uncertainty should widen the envelope"


def test_end_to_end_low_responders_exist():
    """With drug-side uncertainty, some modeled bodies have Emax/ED50 draws
    that leave prolactin below the maturation gate -> zero supply."""
    p = LactationProtocol(galactagogue="metoclopramide",
                          galactagogue_dose_mg_day=10)
    pop = simulate_lactation_population(p, n=200, seed=11,
                                        include_drug_uncertainty=True)
    tail_dt = pop.weeks[1] - pop.weeks[0]
    n_tail = max(1, int(round(8.0 / tail_dt)))
    tail = pop.milk[:, -n_tail:].mean(axis=1)
    assert np.any(tail < 10.0), "low-responding bodies should exist"


def test_population_reproducible_with_seed():
    p = LactationProtocol(galactagogue="domperidone",
                          galactagogue_dose_mg_day=40)
    a = simulate_lactation_population(p, n=30, seed=3)
    b = simulate_lactation_population(p, n=30, seed=3)
    assert np.allclose(a.milk, b.milk)


# --- personalization -----------------------------------------------------------

def _obs(**kw):
    base = dict(weeks_after_drop=4.0, dose_mg_day=10.0,
                prolactin_ng_ml=75.0, milk_ml_day=120.0)
    base.update(kw)
    return Observation(**base)


def test_personalize_runs_and_has_finite_ess():
    res = personalize("metoclopramide", [_obs()], n=150, seed=3)
    assert 0 < res.ess <= 150
    assert np.isfinite(res.ess)
    assert res.weights.sum() == pytest.approx(1.0, rel=1e-9)


def test_personalize_high_prl_at_low_dose_implies_better_responder():
    """Observed prolactin well above the population prediction at 10 mg/d
    should pull the posterior ED50 below the population value (25 mg/d)."""
    res = personalize("metoclopramide", [_obs(prolactin_ng_ml=85.0, milk_ml_day=None)],
                      n=300, seed=3)
    assert res.posterior_param("ed50") < 25.0


def test_personalize_reproducible_with_seed():
    a = personalize("metoclopramide", [_obs()], n=80, seed=9)
    b = personalize("metoclopramide", [_obs()], n=80, seed=9)
    assert np.allclose(a.weights, b.weights)
    assert np.allclose(a.params, b.params)


def test_posterior_predictive_ordered_by_dose():
    res = personalize("metoclopramide", [_obs()], n=120, seed=5)
    proto = LactationProtocol(galactagogue="metoclopramide",
                              galactagogue_dose_mg_day=10)
    txt = posterior_predictive(res, [10, 20], prep_weeks=30,
                               base_protocol=proto)
    # Extract the two median values from the table.
    import re
    meds = [float(m) for m in re.findall(
        r"^\s+\d+\.\d\s*\|\s*\d+\s*\|\s*(\d+)\s*\|", txt, flags=re.M)]
    assert len(meds) == 2
    assert meds[1] >= meds[0], "higher dose -> no lower predicted median"


def test_personalize_requires_observations():
    with pytest.raises(ValueError):
        personalize("metoclopramide", [])
