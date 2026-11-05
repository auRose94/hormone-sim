"""Tests for the low-dose sweep and risk-exposure framing."""

import numpy as np
import pytest

from hormone_sim import (
    dose_sweep_population,
    dose_sweep_table,
    DoseSweepRow,
    exposure_report,
    ExposureProfile,
    PROLACTIN_BASELINE_NG_ML,
)
from hormone_sim.galactagogues import prolactin_response, get_galactagogue


# --- risk-exposure framing ------------------------------------------------------

def test_risk_profile_young_short_course_is_lowest_exposure():
    p = exposure_report(agent="metoclopramide", dose_mg_day=15, weeks=8, age=30)
    assert "low" in p.age_tier
    assert p.cumulative_tier == "low"
    # 8 weeks is a multi-week course: duration risk is accumulating even
    # though it is inside the FDA cap -> 'moderate', honestly.
    assert p.duration_tier.startswith("moderate")


def test_risk_very_short_course_is_low_duration():
    p = exposure_report(agent="metoclopramide", dose_mg_day=10, weeks=3, age=30)
    assert p.duration_tier.startswith("low")


def test_risk_profile_chronic_elderly_is_high():
    p = exposure_report(agent="metoclopramide", dose_mg_day=40, weeks=30, age=72)
    assert p.duration_tier.startswith("HIGH")
    assert p.age_tier.startswith("HIGH")
    assert p.overall_tier == "high"


def test_risk_12_week_boundary_is_flagged():
    """Exactly at the FDA cap the duration tier should be 'moderate', not low."""
    p = exposure_report(agent="metoclopramide", dose_mg_day=15, weeks=12, age=30)
    assert p.duration_tier.startswith("moderate")


def test_risk_beyond_12_weeks_is_high():
    p = exposure_report(agent="metoclopramide", dose_mg_day=15, weeks=16, age=30)
    assert p.duration_tier.startswith("HIGH")


def test_risk_cumulative_is_linear_in_dose_and_time():
    a = exposure_report(agent="metoclopramide", dose_mg_day=10, weeks=10, age=30)
    b = exposure_report(agent="metoclopramide", dose_mg_day=20, weeks=10, age=30)
    assert b.cumulative_mg == pytest.approx(2 * a.cumulative_mg)
    c = exposure_report(agent="metoclopramide", dose_mg_day=10, weeks=20, age=30)
    assert c.cumulative_mg == pytest.approx(2 * a.cumulative_mg)


def test_risk_report_contains_monitoring_and_honesty():
    p = exposure_report(agent="metoclopramide", dose_mg_day=15, weeks=12, age=30)
    txt = p.report()
    assert "AIMS" in txt, "must name the clinical exam (no biomarker exists)"
    assert "not a probability" in txt, "must not pretend to be a risk model"
    assert "post" in txt.lower(), "post-stop monitoring must be mentioned"


def test_risk_concurrent_antipsychotics_flagged():
    p = exposure_report(agent="metoclopramide", dose_mg_day=10, weeks=8,
                        age=30, concurrent_antipsychotics=True)
    assert any("dopamine antagonist" in f for f in p.factors)


def test_risk_unknown_agent_raises():
    with pytest.raises(KeyError):
        exposure_report(agent="fenugreek")  # herbal: not in registry


# --- low-dose sweep ---------------------------------------------------------------

def test_sweep_rows_have_monotone_prolactin():
    rows = dose_sweep_population("metoclopramide", [5, 10, 20, 45], n=30)
    prls = [r.prolactin_ng_ml for r in rows]
    assert all(prls[i] <= prls[i + 1] for i in range(len(prls) - 1))


def test_sweep_prolactin_matches_dose_response():
    rows = dose_sweep_population("metoclopramide", [15], n=10)
    g = get_galactagogue("metoclopramide")
    assert rows[0].prolactin_ng_ml == pytest.approx(
        prolactin_response(g, 15.0), rel=1e-6)


def test_sweep_is_paired_across_doses():
    """Same seed per dose -> the same bodies are evaluated; the population
    medians should move smoothly (no re-shuffling noise between doses)."""
    rows = dose_sweep_population("metoclopramide", [10, 20], n=40, seed=5)
    assert rows[1].milk_p50 >= rows[0].milk_p50


def test_sweep_heterogeneity_shrinks_with_dose():
    """Low doses are bimodal (many bodies at 0, some doing fine); high doses
    compress everyone above the threshold. The compression that matters:
    at the high dose even the 10th-percentile body produces meaningful milk;
    at the low dose it produces nothing."""
    rows = dose_sweep_population("metoclopramide", [5, 30], n=80, seed=3)
    lo, hi = rows[0], rows[-1]
    assert hi.frac_reaching_threshold > lo.frac_reaching_threshold
    assert hi.milk_p10 > 50.0, "high dose: even p10 body reaches meaningful supply"
    assert lo.milk_p10 == 0.0, "low dose: the marginal body produces nothing"


def test_sweep_knee_exists():
    """The model's core low-dose claim: big gains up to ~20 mg/d, then
    diminishing returns. Median milk at 20 should be >= ~85% of the median
    at 45, while exposure (dose) is less than half."""
    rows = dose_sweep_population("metoclopramide", [20, 45], n=80, seed=3)
    m20 = rows[0].milk_p50
    m45 = rows[-1].milk_p50
    assert m20 >= 0.85 * m45


def test_sweep_no_negative_milk():
    rows = dose_sweep_population("metoclopramide", [5, 15, 30], n=30)
    for r in rows:
        assert r.milk_p10 >= 0 and r.milk_p50 >= 0 and r.milk_p90 >= 0


def test_sweep_table_runs():
    rows = dose_sweep_population("metoclopramide", [10, 20], n=20)
    txt = dose_sweep_table(rows)
    assert "heterogeneity" in txt
    assert "metoclopramide" in txt
