"""Tests for estrogenic activity, compounds, and breast development.

These assert the *receptor-level* logic is sound: potency ordering, EE
calculation, inverse-design scaling, and breast development monotonicity.
Clinical accuracy of the calibrated constants is documented as uncertain.
"""

import numpy as np
import pytest

from hormone_sim import (
    PHYTOESTROGENS,
    get_compound,
    hypothetical,
    ee_of,
    required_daily_rate_ug_day,
    TargetProfile,
    simulate_breast_development,
    BreastModel,
)
from hormone_sim.estrogenic_activity import design_summary, format_dose


# --- compounds ---------------------------------------------------------------

def test_known_compounds_present():
    for name in ("genistein", "daidzein", "equol", "coumestrol"):
        c = get_compound(name)
        assert c.rea > 0
        assert c.mw > 0


def test_unknown_compound_raises():
    with pytest.raises(KeyError):
        get_compound("not_a_real_compound")


def test_potency_positive_and_ordered():
    # All potencies positive.
    for c in PHYTOESTROGENS.values():
        assert c.potency > 0
    # More potent REA -> higher potency (same MW/fraction family).
    assert get_compound("coumestrol").potency > get_compound("genistein").potency
    assert get_compound("genistein").potency > get_compound("daidzein").potency


def test_hypothetical_compound():
    c = hypothetical("future", rea=0.5, mw=272.38)
    assert c.name == "future"
    assert c.rea == 0.5
    assert c.potency > 0


def test_potency_scales_with_rea():
    # potency = rea * (MW_e2/MW) * (ff/ff_e2); doubling REA doubles potency.
    a = hypothetical("a", rea=0.01)
    b = hypothetical("b", rea=0.02)
    assert b.potency == pytest.approx(2 * a.potency, rel=1e-9)


# --- estradiol equivalents ---------------------------------------------------

def test_ee_of_straight_line():
    c = get_compound("genistein")
    curve = np.array([100.0, 200.0, 300.0])
    res = ee_of(c, curve)
    expected = curve * c.potency
    assert np.allclose(res.ee_curve_pgml, expected)


def test_ee_of_stronger_compound_higher_ee():
    e2 = hypothetical("e2-ref", rea=1.0, mw=272.38)  # rea=1 -> potency ~1
    gen = get_compound("genistein")
    curve = np.array([100.0])
    ee_e2 = ee_of(e2, curve).ee_curve_pgml[0]
    ee_gen = ee_of(gen, curve).ee_curve_pgml[0]
    assert ee_e2 > ee_gen, "estradiol-equivalent of estradiol > genistein at same pg/mL"


# --- inverse design ----------------------------------------------------------

def test_required_dose_inversely_proportional_to_potency():
    # A weaker compound needs a proportionally larger dose for the same EE.
    strong = hypothetical("strong", rea=0.1)
    weak = hypothetical("weak", rea=0.01)
    d_strong = required_daily_rate_ug_day(strong, 150.0, "transdermal")
    d_weak = required_daily_rate_ug_day(weak, 150.0, "transdermal")
    # weak is 10x weaker -> needs ~10x dose.
    assert d_weak == pytest.approx(10 * d_strong, rel=0.05)


def test_required_dose_positive():
    c = hypothetical("x", rea=0.05)
    d = required_daily_rate_ug_day(c, 150.0, "transdermal")
    assert d > 0


def test_design_summary_runs():
    c = hypothetical("demo", rea=0.01)
    txt = design_summary(c, TargetProfile(ee_target_pgml=150.0))
    assert "demo" in txt
    assert "transdermal" in txt


def test_format_dose_units():
    assert "mg" in format_dose(5000)      # large -> mg
    assert "µg" in format_dose(50)        # small -> µg


# --- breast development ------------------------------------------------------

def test_breast_grows_monotonically():
    ee = np.full(400, 150.0)  # flat EE exposure
    res = simulate_breast_development(ee, weeks=100)
    assert np.all(np.diff(res.volume) >= 0), "volume should be non-decreasing"
    assert res.volume[-1] > res.volume[0], "should grow"


def test_breast_saturates_below_ceiling():
    ee = np.full(1000, 150.0)
    res = simulate_breast_development(ee, weeks=300)
    assert res.volume[-1] <= BreastModel().v_max + 1e-9


def test_breast_higher_ee_faster_growth():
    # 2 years at high EE should reach a higher Tanner stage than low EE.
    hi = simulate_breast_development(np.full(400, 150.0), weeks=100)
    lo = simulate_breast_development(np.full(400, 50.0), weeks=100)
    assert hi.volume[-1] > lo.volume[-1]
    assert hi.final_tanner() >= lo.final_tanner()


def test_breast_tanner_in_range():
    res = simulate_breast_development(np.full(400, 150.0), weeks=100)
    for stage in res.tanner:
        assert 1 <= stage <= 5
