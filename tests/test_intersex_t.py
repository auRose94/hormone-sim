"""Tests for intersex condition packs and trans-masculine personalization."""

import numpy as np
import pytest

from hormone_sim import (
    INTERSEX_PACKS,
    get_pack,
    pack_summary,
    simulate_androgen,
    get_androgen,
    simulate_masculinization,
    AndrogenObservation,
    personalize_androgen,
    posterior_predictive_androgen,
    load_androgen_observations,
    Regimen,
    EstradiolSpec,
    simulate,
)
from hormone_sim.ingest import IngestError


# --- registry -------------------------------------------------------------------

def test_packs_present_and_documented():
    for key in ("cais", "pais", "cah_21ohd", "klinefelter", "turner"):
        p = get_pack(key)
        assert p.citations.strip(), f"{key} missing citations"
        assert p.limits.strip(), f"{key} must state its limits"
        assert p.anchor_quality in ("well-characterized", "moderate",
                                    "framework-only")


def test_unknown_pack_raises():
    with pytest.raises(KeyError):
        get_pack("5alpha_reductase")  # real condition; no DHT layer, not packed


def test_summary_runs():
    txt = pack_summary()
    assert "cais" in txt and "anchor" in txt


# --- pack-applied androgen model ---------------------------------------------------

def test_cah_elevated_endogenous_t():
    """CAH: zero exogenous dose sits at the adrenal-driven baseline (~150),
    not the population ovarian baseline (~30)."""
    a = get_androgen("testosterone_cypionate")
    pack = get_pack("cah_21ohd")
    r = simulate_androgen(a, 0, 7, weeks=6,
                          baseline_t_ngdl=pack.baseline_t_ngdl)
    assert r["t_ngdl"][-1] == pytest.approx(pack.baseline_t_ngdl, rel=0.15)


def test_turner_no_gonadal_output():
    """Turner: zero-dose T sits near zero (streak gonads)."""
    a = get_androgen("testosterone_cypionate")
    pack = get_pack("turner")
    r = simulate_androgen(a, 0, 7, weeks=6,
                          baseline_t_ngdl=pack.baseline_t_ngdl)
    assert r["t_ngdl"][-1] < 15


def test_klinefelter_low_male_baseline():
    """Klinefelter: zero-dose T sits at its documented low-male baseline
    (~250 ng/dL) — above female-typical, below the male reference floor."""
    a = get_androgen("testosterone_cypionate")
    pack = get_pack("klinefelter")
    r = simulate_androgen(a, 0, 7, weeks=6,
                          baseline_t_ngdl=pack.baseline_t_ngdl)
    assert r["t_ngdl"][-1] == pytest.approx(pack.baseline_t_ngdl, rel=0.15)
    assert 100 < r["t_ngdl"][-1] < 300  # low-male band, not female, not full male


def test_pack_masc_passes_through():
    r = simulate_masculinization(
        get_androgen("testosterone_cypionate"), 100, 7, weeks=26,
        baseline_t_ngdl=get_pack("cah_21ohd").baseline_t_ngdl)
    assert r.t_mean_ngdl[-1] > 800  # exogenous dominates; runs fine


# --- pack-applied trans-feminine estradiol model ------------------------------------

def test_cais_zero_dose_sits_at_adrenal_e2():
    """Post-gonadectomy CAIS: zero exogenous E2 settles at adrenal-only
    levels (~8 pg/mL), not the AMAB population baseline (15)."""
    pack = get_pack("cais")
    r = Regimen(name="cais0", population="trans_feminine",
                estradiol=EstradiolSpec(route="transdermal", daily_rate_ug_day=0))
    sim = simulate(r, weeks=8, baseline_e2_pgml=pack.tf_baseline_e2_pgml)
    assert sim.e2_pgml[-1] == pytest.approx(pack.tf_baseline_e2_pgml, rel=0.2)


def test_cais_responds_normally_to_e2():
    """CAIS: ER intact — the estradiol dose-response is standard."""
    pack = get_pack("cais")
    r = Regimen(name="cais50", population="trans_feminine",
                estradiol=EstradiolSpec(route="transdermal", daily_rate_ug_day=50))
    sim = simulate(r, weeks=12, baseline_e2_pgml=pack.tf_baseline_e2_pgml)
    assert 100 <= sim.steady_state_e2_pgml() <= 220  # same target band


def test_turner_trans_fem_baseline():
    pack = get_pack("turner")
    r = Regimen(name="turner0", population="trans_feminine",
                estradiol=EstradiolSpec(route="transdermal", daily_rate_ug_day=0))
    sim = simulate(r, weeks=8, baseline_e2_pgml=pack.tf_baseline_e2_pgml)
    assert sim.e2_pgml[-1] == pytest.approx(6.0, rel=0.2)


# --- T-log ingest --------------------------------------------------------------------

@pytest.fixture
def t_log(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text(
        "date,start_date,weeks_on_t,dose_mg,interval_days,"
        "testosterone_ngdl,hct_pct,notes\n"
        "2026-08-07,2026-08-07,,100,7,410,41.2,\n"
        ",2026-08-07,6,100,7,880,,\n"
        ",2026-08-07,12,100,7,,45.8,cbc only\n"
    )
    return p


def test_t_ingest_mixed_axes(t_log):
    obs, skipped = load_androgen_observations(t_log)
    assert len(obs) == 3 and skipped == 0
    assert obs[0].weeks_on_t == pytest.approx(0.0)
    assert obs[2].hct_pct == 45.8 and obs[2].testosterone_ngdl is None


def test_t_ingest_fails_closed(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("weeks_on_t,dose_mg,testosterone_ngdl\n"
                 "4,,900\n")  # missing dose, no default
    with pytest.raises(IngestError):
        load_androgen_observations(p)


def test_t_ingest_default_interval(tmp_path):
    p = tmp_path / "noint.csv"
    p.write_text("weeks_on_t,dose_mg,testosterone_ngdl\n6,100,880\n")
    obs, _ = load_androgen_observations(p, default_interval_days=7.0)
    assert obs[0].interval_days == 7.0


# --- trans-masculine personalization ---------------------------------------------------

def _obs(**kw):
    base = dict(weeks_on_t=6.0, dose_mg=100.0, interval_days=7.0,
                testosterone_ngdl=880.0, hct_pct=None)
    base.update(kw)
    return AndrogenObservation(**base)


def test_personalize_androgen_runs_and_reproducible():
    a = personalize_androgen("testosterone_cypionate", [_obs()], n=60, seed=3)
    b = personalize_androgen("testosterone_cypionate", [_obs()], n=60, seed=3)
    assert 0 < a.ess <= 60
    assert np.allclose(a.weights, b.weights)
    assert a.weights.sum() == pytest.approx(1.0, rel=1e-9)


def test_personalize_androgen_requires_observations():
    with pytest.raises(ValueError):
        personalize_androgen("testosterone_cypionate", [])


def test_high_t_observation_pulls_vd_down():
    """Measured T well above prior prediction at the same dose implies a
    smaller volume (less dilution). Directional sanity of the inference."""
    res = personalize_androgen(
        "testosterone_cypionate",
        [_obs(testosterone_ngdl=1400.0, hct_pct=None)], n=200, seed=3)
    assert res.posterior_param("vd") < 2.0


def test_predictive_higher_dose_higher_trough():
    res = personalize_androgen("testosterone_cypionate", [_obs()], n=60, seed=3)
    txt = posterior_predictive_androgen(res, "testosterone_cypionate",
                                        [(100, 7), (150, 7)])
    import re
    troughs = [float(m) for m in re.findall(
        r"^\s+\d+\.\d\s*\|\s+\d+\.\d\s*\|\s*(\d+)\s*\|", txt, flags=re.M)]
    assert len(troughs) == 2
    assert troughs[1] > troughs[0]


def test_pack_applied_personalization():
    """CAH participant: the pack's adrenal baseline is the prior's center."""
    res = personalize_androgen(
        "testosterone_cypionate", [_obs()], n=60, seed=3,
        baseline_t_ngdl=get_pack("cah_21ohd").baseline_t_ngdl)
    assert res.posterior_param("base_t") > 60  # centered far above ovarian
