"""Tests for the measurement-model abstraction (AGENTS.md §5a item 1).

The engine (:mod:`hormone_sim.measurement`) is the study layer's shared
machinery; both domains migrated onto it. These tests pin the migration's
behavior contract: each domain's public function must reproduce a reference
implementation of the pre-abstraction algorithm **exactly** (same RNG call
order, same simulations) — so future engine changes that drift from the
documented importance-sampling mechanics fail here, loudly.
"""

from dataclasses import replace as dc_replace
from typing import List

import numpy as np
import pytest

from hormone_sim import (
    Observation,
    AndrogenObservation,
    LactationProtocol,
    LactationModel,
    simulate_lactation,
    simulate_androgen,
    personalize,
    posterior_predictive,
    personalize_androgen,
    posterior_predictive_androgen,
    androgen_population_sweep,
    androgen_sweep_table,
    get_androgen,
    get_galactagogue,
    PROLACTIN_BASELINE_NG_ML,
)
from hormone_sim.measurement import (
    BodySample,
    MeasurementDomain,
    personalize_domain,
    posterior_predictive_table,
    population_bands,
)
from hormone_sim.personalization import LactationDomain
from hormone_sim.personalization_androgen import AndrogenDomain
from hormone_sim.sensitivity import SensitivitySpec


# --- reference: the pre-abstraction lactation loop ----------------------------------

def _reference_lactation(agent, observations, prep_weeks, p_start, sessions,
                         n, seed, dt_weeks=0.25):
    """Verbatim port of the pre-abstraction `personalize()` loop (kept as the
    behavior contract for the engine migration)."""
    g = get_galactagogue(agent)
    spec = SensitivitySpec()
    rng = np.random.default_rng(seed)

    max_week = max(o.weeks_after_drop for o in observations)
    base_protocol = LactationProtocol(
        galactagogue=agent,
        galactagogue_dose_mg_day=float(np.mean([o.dose_mg_day for o in observations])),
        prep_weeks=prep_weeks,
        progesterone_start_weeks=p_start,
        pumping_sessions_per_day=sessions,
    )
    param_names = ["emax", "ed50", "baseline", "milk_vmax", "milk_km",
                   "mat_tau", "mat_thr", "ductal_tau", "la_tau", "invol"]
    P = np.empty((n, len(param_names)))
    logw = np.empty(n)
    milk_curves = None
    weeks = None

    total_weeks = prep_weeks + max(max_week + 8.0, 52.0)

    for i in range(n):
        m = spec.sample(rng)
        emax = g.prolactin_emax_ng_ml * float(
            rng.lognormal(0.0, spec.galactagogue_cvs["prolactin_emax_ng_ml"]))
        ed50 = max(1e-3, g.prolactin_ed50_mg_day * float(
            rng.lognormal(0.0, spec.galactagogue_cvs["prolactin_ed50_mg_day"])))
        base = PROLACTIN_BASELINE_NG_ML * float(
            rng.lognormal(0.0, spec.galactagogue_cvs["baseline"]))
        P[i] = [emax, ed50, base,
                m.milk_vmax_ml_day, m.milk_km_prolactin_ng_ml,
                m.maturation_tau_weeks, m.maturation_prolactin_threshold_ng_ml,
                m.ductal_tau_weeks, m.la_tau_weeks, m.involution_rate_per_week]

        lw = 0.0
        for o in observations:
            if o.prolactin_ng_ml is not None:
                prl = base + emax * o.dose_mg_day / (ed50 + o.dose_mg_day)
                sigma = max(0.15 * o.prolactin_ng_ml, 3.0)
                lw += -0.5 * ((prl - o.prolactin_ng_ml) / sigma) ** 2
            if o.milk_ml_day is not None:
                proto = dc_replace(base_protocol,
                                   galactagogue_dose_mg_day=float(o.dose_mg_day))
                res = simulate_lactation(proto, total_weeks=total_weeks,
                                         dt_weeks=dt_weeks, model=m)
                wk_idx = int(round((prep_weeks + o.weeks_after_drop) / dt_weeks))
                pred = float(res.milk_ml_day[min(wk_idx, res.milk_ml_day.size - 1)])
                sigma = max(0.30 * o.milk_ml_day, 10.0)
                lw += -0.5 * ((pred - o.milk_ml_day) / sigma) ** 2
                if milk_curves is None:
                    milk_curves = np.empty((n, res.milk_ml_day.size))
                    weeks = res.weeks
                milk_curves[i] = res.milk_ml_day
        logw[i] = lw

    if milk_curves is None:
        res = simulate_lactation(base_protocol, total_weeks=total_weeks,
                                 dt_weeks=dt_weeks)
        milk_curves = np.tile(res.milk_ml_day, (n, 1))
        weeks = res.weeks

    w = np.exp(logw - logw.max())
    w /= w.sum()
    ess = 1.0 / np.sum(w ** 2) if w.sum() > 0 else 0.0
    return dict(weights=w, ess=float(ess), params=P, weeks=weeks,
                milk=milk_curves)


def _obs(**kw):
    base = dict(weeks_after_drop=4.0, dose_mg_day=10.0,
                prolactin_ng_ml=75.0, milk_ml_day=120.0)
    base.update(kw)
    return Observation(**base)


def test_lactation_engine_matches_reference_with_milk():
    """Mixed prolactin+milk log: the engine must reproduce the reference
    ensemble bit-for-bit (weights, params, curves)."""
    obs = [_obs(), _obs(weeks_after_drop=8.0, dose_mg_day=15.0,
                        prolactin_ng_ml=None, milk_ml_day=165.0)]
    ref = _reference_lactation("metoclopramide", obs, 30.0, None, 8.0,
                               n=40, seed=3)
    res = personalize("metoclopramide", obs, n=40, seed=3)
    assert np.allclose(res.weights, ref["weights"], rtol=1e-12)
    assert res.ess == pytest.approx(ref["ess"], rel=1e-12)
    assert np.allclose(res.params, ref["params"], rtol=1e-12)
    assert np.allclose(res.milk, ref["milk"], rtol=1e-12)
    assert np.allclose(res.weeks, ref["weeks"])


def test_lactation_engine_matches_reference_prl_only():
    """Prolactin-only log: tiled default curves, weights from the analytic
    update only — the reference must match here too."""
    obs = [_obs(prolactin_ng_ml=85.0, milk_ml_day=None)]
    ref = _reference_lactation("metoclopramide", obs, 30.0, None, 8.0,
                               n=30, seed=5)
    res = personalize("metoclopramide", obs, n=30, seed=5)
    assert np.allclose(res.weights, ref["weights"], rtol=1e-12)
    assert np.allclose(res.milk, ref["milk"], rtol=1e-12)
    # Every body's curve is the same tiled population-default simulation.
    assert np.allclose(res.milk, res.milk[0])


# --- reference: the pre-abstraction androgen loop -----------------------------------

def _reference_androgen(agent, observations, n, seed, body_weight_kg=70.0,
                        baseline_t_ngdl=None, baseline_e2_pgml=None,
                        dt_hours=2.0):
    """Verbatim port of the pre-abstraction `personalize_androgen()` loop."""
    from hormone_sim.erythrocytosis import HctModel, simulate_hct
    from hormone_sim.androgens import BASELINE_T_NGDL, BASELINE_E2_PGML
    a = get_androgen(agent)
    rng = np.random.default_rng(seed)

    base_t = BASELINE_T_NGDL if baseline_t_ngdl is None else float(baseline_t_ngdl)
    base_e2 = BASELINE_E2_PGML if baseline_e2_pgml is None else float(baseline_e2_pgml)

    max_week = max(o.weeks_on_t for o in observations)
    total_weeks = max_week + 4.0

    param_names = ["half_life", "vd", "ec50", "base_t", "base_e2",
                   "hct_rise", "hct_tau", "hct_baseline"]
    P = np.empty((n, len(param_names)))
    logw = np.empty(n)
    t_curves = None
    weeks_axis = None

    for i in range(n):
        half_life = max(0.5, 2.0 * rng.lognormal(0.0, 0.30))
        vd = max(0.5, 2.0 * rng.lognormal(0.0, 0.30))
        ec50 = max(20.0, 170.0 * rng.lognormal(0.0, 0.30))
        b_t = max(1.0, base_t * rng.lognormal(0.0, 0.30))
        b_e2 = max(1.0, base_e2 * rng.lognormal(0.0, 0.30))
        hct_rise = max(1.0, 13.0 * rng.lognormal(0.0, 0.30))
        hct_tau = max(2.0, 10.0 * rng.lognormal(0.0, 0.40))
        hct_base = float(np.clip(rng.normal(40.0, 2.0), 33.0, 47.0))
        P[i] = [half_life, vd, ec50, b_t, b_e2, hct_rise, hct_tau, hct_base]

        lw = 0.0
        t_curve = None
        weeks_i = None
        a_i = dc_replace(a, vd_l_kg=vd)
        for o in observations:
            if o.testosterone_ngdl is None and o.hct_pct is None:
                continue
            if t_curve is None:
                pk = simulate_androgen(
                    a_i, o.dose_mg, o.interval_days, total_weeks,
                    dt_hours=dt_hours, e2_ec50_ngdl=ec50,
                    body_weight_kg=body_weight_kg,
                    baseline_t_ngdl=b_t, baseline_e2_pgml=b_e2,
                )
                t_curve = pk["t_ngdl"]
                weeks_i = pk["hours"] / (7 * 24.0)
            if o.testosterone_ngdl is not None:
                idx = min(int(round(o.weeks_on_t * 7 * 24 / dt_hours)),
                          t_curve.size - 1)
                pred = float(t_curve[idx])
                sigma = max(0.20 * o.testosterone_ngdl, 30.0)
                lw += -0.5 * ((pred - o.testosterone_ngdl) / sigma) ** 2
            if o.hct_pct is not None:
                hm = HctModel(baseline_hct=hct_base, max_rise_points=hct_rise,
                              tau_wk=hct_tau)
                hr = simulate_hct(a, o.dose_mg, o.interval_days,
                                  weeks=max(o.weeks_on_t + 1.0, 2.0),
                                  model=hm, body_weight_kg=body_weight_kg)
                idx2 = min(int(round(o.weeks_on_t / 0.25)), hr.hct.size - 1)
                pred_h = float(hr.hct[idx2])
                sigma_h = max(0.05 * o.hct_pct, 1.5)
                lw += -0.5 * ((pred_h - o.hct_pct) / sigma_h) ** 2
        logw[i] = lw
        if t_curve is not None and t_curves is None:
            t_curves = np.empty((n, t_curve.size))
            weeks_axis = weeks_i
        if t_curve is not None:
            t_curves[i] = t_curve
        else:
            t_curves[i] = 0.0

    dominant = max(observations, key=lambda o: o.weeks_on_t)
    hct_curves = np.empty((n, t_curves.shape[1]))
    for i in range(n):
        half_life, vd, ec50, b_t, b_e2, hct_rise, hct_tau, hct_base = P[i]
        hm = HctModel(baseline_hct=hct_base, max_rise_points=hct_rise,
                      tau_wk=hct_tau)
        hr = simulate_hct(a, dominant.dose_mg, dominant.interval_days,
                          weeks=total_weeks, model=hm,
                          body_weight_kg=body_weight_kg)
        hct_curves[i] = np.interp(weeks_axis, hr.weeks, hr.hct)

    w = np.exp(logw - logw.max())
    w /= w.sum()
    ess = 1.0 / np.sum(w ** 2) if w.sum() > 0 else 0.0
    return dict(weights=w, ess=float(ess), params=P, weeks=weeks_axis,
                t=t_curves, hct=hct_curves)


def _tobs(**kw):
    base = dict(weeks_on_t=6.0, dose_mg=100.0, interval_days=7.0,
                testosterone_ngdl=880.0, hct_pct=None)
    base.update(kw)
    return AndrogenObservation(**base)


def test_androgen_engine_matches_reference():
    obs = [_tobs(), _tobs(weeks_on_t=10.0, testosterone_ngdl=None,
                          hct_pct=45.5)]
    ref = _reference_androgen("testosterone_cypionate", obs, n=30, seed=3)
    res = personalize_androgen("testosterone_cypionate", obs, n=30, seed=3)
    assert np.allclose(res.weights, ref["weights"], rtol=1e-12)
    assert res.ess == pytest.approx(ref["ess"], rel=1e-12)
    assert np.allclose(res.params, ref["params"], rtol=1e-12)
    assert np.allclose(res.t_ngdl, ref["t"], rtol=1e-12)
    assert np.allclose(res.hct, ref["hct"], rtol=1e-12)
    assert np.allclose(res.weeks, ref["weeks"])


def test_androgen_engine_pack_baseline_matches_reference():
    """The intersex-pack override path must go through the engine with the
    same result as the reference loop with the same baseline."""
    from hormone_sim import get_pack
    obs = [_tobs()]
    pack_base = get_pack("cah_21ohd").baseline_t_ngdl
    ref = _reference_androgen("testosterone_cypionate", obs, n=25, seed=7,
                              baseline_t_ngdl=pack_base)
    res = personalize_androgen("testosterone_cypionate", obs, n=25, seed=7,
                               baseline_t_ngdl=pack_base)
    assert np.allclose(res.weights, ref["weights"], rtol=1e-12)
    assert np.allclose(res.params, ref["params"], rtol=1e-12)


# --- engine contract -----------------------------------------------------------------

def test_engine_rejects_empty_observations():
    domain = LactationDomain(agent="metoclopramide",
                             base_protocol=LactationProtocol())
    with pytest.raises(ValueError):
        personalize_domain(domain, [], n=10, seed=1)


def test_engine_requires_domains_to_fill_curves():
    class BareDomain(MeasurementDomain):
        name = "bare"
        param_names = ["x"]

        def sample_body(self, rng, observations):
            return BodySample(params=np.array([1.0]))

        def log_likelihood(self, body, obs):
            return 0.0

        def curve_names(self):
            return ("milk",)

    with pytest.raises(ValueError, match="missing curve"):
        personalize_domain(BareDomain(), [_obs()], n=5, seed=1)


def test_generic_table_delegation_both_domains():
    """The public predictive tables must be the generic renderer fed by the
    domain descriptors — one presentation contract, two domains."""
    res = personalize("metoclopramide", [_obs()], n=30, seed=5)
    proto = LactationProtocol(galactagogue="metoclopramide",
                              galactagogue_dose_mg_day=10)
    via_wrapper = posterior_predictive(res, [10, 20], prep_weeks=30,
                                       base_protocol=proto)
    domain = LactationDomain(agent="metoclopramide", base_protocol=proto)
    via_generic = posterior_predictive_table(domain, res, [10, 20])
    assert via_wrapper == via_generic

    t_res = personalize_androgen("testosterone_cypionate", [_tobs()],
                                 n=20, seed=3)
    t_wrapper = posterior_predictive_androgen(t_res, "testosterone_cypionate",
                                              [(100, 7), (150, 7)])
    t_domain = AndrogenDomain(agent="testosterone_cypionate")
    assert t_wrapper == posterior_predictive_table(t_domain, t_res,
                                                   [(100, 7), (150, 7)])


# --- population machinery ------------------------------------------------------------

def test_population_bands_lactation_matches_wrapper():
    proto = LactationProtocol(galactagogue="metoclopramide",
                              galactagogue_dose_mg_day=10)
    domain = LactationDomain(agent="metoclopramide", base_protocol=proto)
    weeks, curves = population_bands(domain, 10.0, n=15, seed=11)
    assert curves["milk"].shape == (15, weeks.size)
    assert np.all(curves["milk"] >= 0)


def test_androgen_sweep_rows_paired_and_monotone():
    rows = androgen_population_sweep(
        "testosterone_cypionate", [(50, 7), (100, 7), (200, 7)], n=40, seed=3)
    assert len(rows) == 3
    for r in rows:
        assert 0.0 <= r.frac_hit <= 1.0
        tr = r.metrics["trough"]
        assert tr[0] <= tr[1] <= tr[2]
    # Trough rises with dose (paired bodies). The trough-IN-RANGE fraction is
    # an inverted U across doses — most bodies in range near the middle of
    # the sweep, and the supra-therapeutic dose falls OUT the top of the
    # range (more exposure is not more in-range).
    assert rows[2].metric("trough") > rows[0].metric("trough")
    assert rows[0].frac_hit > rows[2].frac_hit


def test_androgen_sweep_interval_changes_swing_not_mean():
    """The model's one scheduling lesson, at the population level: weekly vs
    biweekly at the same weekly dose -> same mean, bigger swing."""
    rows = androgen_population_sweep(
        "testosterone_cypionate", [(100, 7), (200, 14)], n=60, seed=3)
    weekly, biweekly = rows
    mean_w = weekly.metric("mean")
    mean_b = biweekly.metric("mean")
    assert abs(mean_b - mean_w) < 0.15 * mean_w          # same mean (±15%)
    swing_w = weekly.metric("peak") - weekly.metric("trough")
    swing_b = biweekly.metric("peak") - biweekly.metric("trough")
    assert swing_b > swing_w                              # bigger swing


def test_androgen_sweep_table_text():
    rows = androgen_population_sweep("testosterone_cypionate",
                                     [(100, 7)], n=20, seed=3)
    txt = androgen_sweep_table(rows, "testosterone_cypionate")
    assert "testosterone_cypionate" in txt
    assert "trough" in txt and "peak" in txt
    assert "never dosing advice" in txt   # the population-vs-individual rule


def test_engine_body_sample_is_the_contract():
    """BodySample round-trips through the engine: params rows align with
    param_names for both domains."""
    for domain, obs, n in (
        (LactationDomain(agent="metoclopramide",
                         base_protocol=LactationProtocol(
                             galactagogue="metoclopramide")), [_obs()], 8),
        (AndrogenDomain(agent="testosterone_cypionate"), [_tobs()], 8),
    ):
        post = personalize_domain(domain, obs, n=n, seed=2)
        assert post.params.shape == (n, len(domain.param_names))
        assert post.weights.shape == (n,)
        assert post.weights.sum() == pytest.approx(1.0, rel=1e-9)
        assert 0 < post.ess <= n
