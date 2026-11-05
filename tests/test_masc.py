"""Tests for the trans-masculine domain: androgen PK, masculinization PD,
erythrocytosis framing."""

import numpy as np
import pytest

from hormone_sim import (
    ANDROGENS,
    get_androgen,
    simulate_androgen,
    ugl_to_ngdl,
    ngdl_to_ugl,
    MasculinizationModel,
    simulate_masculinization,
    simulate_masc_population,
    HctModel,
    simulate_hct,
    hct_risk_report,
)


# --- units ----------------------------------------------------------------------

def test_unit_conversions():
    assert ngdl_to_ugl(500.0) == pytest.approx(5.0)
    assert ugl_to_ngdl(5.0) == pytest.approx(500.0)
    assert ugl_to_ngdl(ngdl_to_ugl(750.0)) == pytest.approx(750.0)


# --- androgen registry ------------------------------------------------------------

def test_registered_agents_present():
    for name in ("testosterone_cypionate", "testosterone_enanthate",
                 "testosterone_gel", "testosterone_undecanoate_im"):
        a = get_androgen(name)
        assert a.half_life_h > 0
        assert a.vd_l_kg > 0


def test_oral_undecanoate_deliberately_absent():
    with pytest.raises(KeyError):
        get_androgen("testosterone_undecanoate_oral")


def test_every_agent_has_citation():
    for a in ANDROGENS.values():
        assert a.citation.strip()


def test_depot_rate_positive_for_esters():
    assert get_androgen("testosterone_cypionate").depot_rate_h() > 0
    assert get_androgen("testosterone_gel").depot_rate_h() == 0.0


# --- T PK ---------------------------------------------------------------------------

def _run(name, dose, interval, weeks=12):
    return simulate_androgen(get_androgen(name), dose, interval, weeks=weeks)


def test_cypionate_biweekly_peak_trough_pattern():
    """200 mg/2wk: peak well above the male ref range, trough still male-ish.
    This is the clinically-discussed pattern the model must reproduce."""
    r = _run("testosterone_cypionate", 200, 14, weeks=12)
    ss = r["t_ngdl"][r["t_ngdl"].size // 2:]
    assert ss.max() > 1200, "biweekly peak should ride high"
    assert 300 < ss.min() < 900, "trough should stay in/near male range"
    assert ss.mean() > 800


def test_same_weekly_dose_same_mean_different_swing():
    """100 mg weekly vs 200 mg biweekly: same average, tighter swing weekly."""
    wk = _run("testosterone_cypionate", 100, 7, weeks=12)
    bi = _run("testosterone_cypionate", 200, 14, weeks=12)
    w_ss = wk["t_ngdl"][wk["t_ngdl"].size // 2:]
    b_ss = bi["t_ngdl"][bi["t_ngdl"].size // 2:]
    assert abs(w_ss.mean() - b_ss.mean()) < 0.25 * b_ss.mean()
    swing_w = w_ss.max() - w_ss.min()
    swing_b = b_ss.max() - b_ss.min()
    assert swing_w < swing_b, "weekly should swing less than biweekly"


def test_gel_is_flat():
    r = _run("testosterone_gel", 50, 1, weeks=10)
    ss = r["t_ngdl"][r["t_ngdl"].size // 2:]
    assert ss.max() - ss.min() < 0.05 * ss.mean()


def test_gel_lower_peak_than_biweekly_at_typical_doses():
    gel = _run("testosterone_gel", 50, 1, weeks=10)
    bi = _run("testosterone_cypionate", 200, 14, weeks=10)
    assert gel["t_ngdl"].max() < bi["t_ngdl"][bi["t_ngdl"].size // 2:].max()


def test_exogenous_t_suppresses_ovarian_output():
    """The reversed-direction HPG feedback: ovarian fraction falls."""
    r = _run("testosterone_cypionate", 100, 7, weeks=10)
    assert r["ovarian_frac"][-1] < 0.25
    assert r["ovarian_frac"][0] > 0.9  # starts intact


def test_e2_falls_into_moderate_range_on_t():
    """Aromatized E2 on T: tens of pg/mL, far below ovarian baseline peaks."""
    r = _run("testosterone_cypionate", 100, 7, weeks=10)
    e2 = r["e2_pgml"][-1]
    assert 10 < e2 < 90


def test_more_dose_higher_t():
    lo = _run("testosterone_cypionate", 50, 7, weeks=10)
    hi = _run("testosterone_cypionate", 150, 7, weeks=10)
    assert hi["t_ngdl"][-1] > lo["t_ngdl"][-1]


# --- masculinization PD -----------------------------------------------------------

def test_cycle_cessation_anchor():
    """Median cessation ~5 months (22 wk) at a typical regimen (the anchor)."""
    r = simulate_masculinization(get_androgen("testosterone_cypionate"),
                                 100, 7, weeks=104)
    wk = r.weeks_to_frac(0.5)
    assert wk is not None
    assert 12 <= wk <= 40, f"anchor drifted: 50% at {wk} wk"


def test_cessation_higher_with_more_exposure():
    lo = simulate_masculinization(get_androgen("testosterone_gel"), 25, 1, weeks=104)
    hi = simulate_masculinization(get_androgen("testosterone_cypionate"), 150, 7, weeks=104)
    assert hi.cycle_ceased_prob[-1] >= lo.cycle_ceased_prob[-1]


def test_outcomes_monotonic_and_ordered():
    """Bottom growth fastest, hair slowest — the clinical ordering."""
    r = simulate_masculinization(get_androgen("testosterone_cypionate"),
                                 100, 7, weeks=104)
    for arr in (r.voice, r.hair, r.body_comp, r.bottom):
        assert np.all(np.diff(arr) >= 0)
    i = r.weeks.size - 1
    assert r.bottom[i] > r.voice[i] > r.hair[i]


def test_report_weeks_are_weeks_not_indices():
    """Regression: the report must index by week, not by 0.25-wk steps."""
    r = simulate_masculinization(get_androgen("testosterone_cypionate"),
                                 100, 7, weeks=104)
    dt = r.weeks[1] - r.weeks[0]
    manual = r.voice[int(round(104 / dt))]
    assert f"{manual:>6.2f}" in r.report(), "wk104 column must match the curve"


def test_population_bands_genetics_spread():
    pop = simulate_masc_population(get_androgen("testosterone_cypionate"),
                                   100, 7, weeks=104, n=60, seed=3)
    lo, med, hi = pop.outcome_band("voice")
    assert lo <= med <= hi
    assert hi > lo, "genetics spread should be real"
    assert np.allclose(pop.cycle_ceased_prob[:, -1].clip(0, 1),
                       pop.cycle_ceased_prob[:, -1])


def test_population_reproducible():
    a = simulate_masc_population(get_androgen("testosterone_cypionate"),
                                 100, 7, weeks=52, n=30, seed=5)
    b = simulate_masc_population(get_androgen("testosterone_cypionate"),
                                 100, 7, weeks=52, n=30, seed=5)
    assert np.allclose(a.voice, b.voice)


# --- erythrocytosis -----------------------------------------------------------------

def test_hct_rises_and_plateaus():
    r = simulate_hct(get_androgen("testosterone_cypionate"), 100, 7, weeks=52)
    assert np.all(np.diff(r.hct) >= -1e-9)
    assert 44 < r.steady_state_hct() < 54


def test_hct_route_dependence():
    """Injectable ester drives a higher plateau than gel — the documented
    route dependence, from peak exposure."""
    inj = simulate_hct(get_androgen("testosterone_cypionate"), 100, 7, weeks=52)
    gel = simulate_hct(get_androgen("testosterone_gel"), 50, 1, weeks=52)
    assert inj.steady_state_hct() > gel.steady_state_hct()


def test_hct_lower_at_lower_dose():
    hi = simulate_hct(get_androgen("testosterone_cypionate"), 200, 14, weeks=52)
    lo = simulate_hct(get_androgen("testosterone_cypionate"), 50, 14, weeks=52)
    assert hi.steady_state_hct() > lo.steady_state_hct()


def test_hct_report_names_the_asymmetry():
    txt = hct_risk_report(get_androgen("testosterone_cypionate"), 100, 7)
    assert "CBC" in txt, "the biomarker point is the whole framing"
    assert "reversible" in txt.lower()
    assert "not a dosing decision" in txt


def test_hct_report_route_note():
    txt = hct_risk_report(get_androgen("testosterone_gel"), 50, 1)
    assert "gel" in txt.lower()
