"""Tests for the study-tool layer: CSV ingest, AIMS tracking, protocol."""

import numpy as np
import pytest

from hormone_sim import (
    load_observations,
    personalize_from_csv,
    IngestError,
    AIMSVisit,
    AIMSRecord,
    parse_visit,
    AIMS_ITEMS,
    generate_protocol,
    ProtocolConfig,
    generate_androgen_protocol,
    AndrogenProtocolConfig,
    Observation,
    personalize,
)


# --- CSV ingest ------------------------------------------------------------------

@pytest.fixture
def log_csv(tmp_path):
    p = tmp_path / "log.csv"
    p.write_text(
        "date,drop_date,weeks_post_drop,dose_mg_day,prolactin_ng_ml,milk_ml_day,notes\n"
        "2026-09-17,2026-09-17,,10,58,,\n"
        ",2026-09-17,2,10,72,95,\n"
        ",2026-09-17,6,10,78,140,\n"
        ",2026-09-17,8,15,,165,\n"
    )
    return p


def test_ingest_reads_mixed_time_axes(log_csv):
    obs, skipped = load_observations(log_csv)
    assert len(obs) == 4 and skipped == 0
    assert obs[0].weeks_after_drop == pytest.approx(0.0)   # from dates
    assert obs[1].weeks_after_drop == pytest.approx(2.0)   # explicit weeks
    assert obs[3].dose_mg_day == 15.0                      # per-row dose


def test_ingest_oz_conversion(tmp_path):
    p = tmp_path / "oz.csv"
    p.write_text("weeks_post_drop,dose_mg_day,milk_ml_day\n"
                 "2,10,3.3814\n")  # ~100 mL
    obs, _ = load_observations(p, milk_unit="oz")
    assert obs[0].milk_ml_day == pytest.approx(100.0, rel=0.01)


def test_ingest_skips_unmeasured_rows(tmp_path):
    p = tmp_path / "gap.csv"
    p.write_text("weeks_post_drop,dose_mg_day,prolactin_ng_ml,milk_ml_day\n"
                 "1,10,60,\n"
                 "2,10,,\n"          # diary-only row -> skipped
                 "3,10,70,\n")
    obs, skipped = load_observations(p)
    assert len(obs) == 2 and skipped == 1


def test_ingest_fails_closed_listing_all_problems(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("weeks_post_drop,dose_mg_day,milk_ml_day\n"
                 "two,10,100\n"      # bad week
                 "3,,100\n"          # missing dose
                 "-1,10,100\n")      # negative week
    with pytest.raises(IngestError) as e:
        load_observations(p)
    msg = str(e.value)
    assert "row 2" in msg and "row 3" in msg and "row 4" in msg


def test_ingest_missing_file():
    with pytest.raises(FileNotFoundError):
        load_observations("/nonexistent/log.csv")


def test_ingest_dose_required(tmp_path):
    p = tmp_path / "nodose.csv"
    p.write_text("weeks_post_drop,milk_ml_day\n2,100\n")
    with pytest.raises(IngestError):
        load_observations(p)
    obs, _ = load_observations(p, default_dose_mg_day=10.0)
    assert obs[0].dose_mg_day == 10.0


def test_personalize_from_csv_matches_direct(log_csv):
    """The CSV path must be equivalent to hand-building the Observations."""
    obs, _ = load_observations(log_csv)
    a = personalize_from_csv("metoclopramide", log_csv, n=50, seed=1)
    b = personalize("metoclopramide", obs, n=50, seed=1)
    assert np.allclose(a.weights, b.weights)


# --- AIMS -------------------------------------------------------------------------

def test_aims_uniform_visit_and_total():
    v = AIMSVisit.uniform(0.0, 0)
    assert v.total() == 0
    v2 = AIMSVisit.uniform(4.0, 1)
    assert v2.total() == 10


def test_aims_rejects_incomplete_exam():
    with pytest.raises(ValueError):
        AIMSVisit(week=4.0, scores={"lips_perioral": 2})


def test_aims_rejects_out_of_range():
    with pytest.raises(ValueError):
        AIMSVisit.uniform(4.0, 5)


def test_aims_flags_item_ge3():
    r = AIMSRecord(baseline=AIMSVisit.uniform(0.0, 0))
    flags = r.add_visit(parse_visit("8:lips_perioral=3"))
    assert any(">=3" in f for f in flags)


def test_aims_flags_rise_above_baseline():
    r = AIMSRecord(baseline=AIMSVisit.uniform(0.0, 1))  # non-zero baseline
    flags = r.add_visit(parse_visit("8:lips_perioral=3"))  # +2 vs baseline
    assert any("above baseline" in f for f in flags)


def test_aims_no_flags_on_stable_record():
    r = AIMSRecord(baseline=AIMSVisit.uniform(0.0, 0))
    flags = r.add_visit(parse_visit("4:all=0"))
    assert flags == []


def test_aims_parse_zero_fill():
    v = parse_visit("8:lips_perioral=2,tongue=1")
    assert len(v.scores) == len(AIMS_ITEMS)
    assert v.scores["jaw"] == 0


def test_aims_parse_all_zero_is_not_falsy():
    """Regression: 'all=0' must build a complete exam, not skip uniform."""
    v = parse_visit("0:all=0")
    assert v.total() == 0 and len(v.scores) == len(AIMS_ITEMS)


def test_aims_post_stop_visit_recorded():
    r = AIMSRecord(baseline=AIMSVisit.uniform(0.0, 0))
    r.add_visit(parse_visit("13:lips_perioral=1,after_stop"))
    assert r.visits[0].after_stop is True
    assert "(post-stop)" in r.report()


def test_aims_nonzero_baseline_note():
    r = AIMSRecord(baseline=AIMSVisit.uniform(0.0, 1))
    assert "Non-zero baseline" in r.notes


# --- protocol skeleton --------------------------------------------------------------

def test_protocol_contains_required_sections():
    doc = generate_protocol(ProtocolConfig(site="Test Clinic", pi="Dr. X"))
    for section in ("Background", "Specific aims", "Design", "Population",
                    "Interventions", "Expected values", "Safety and stopping",
                    "Ethics and consent", "Statistical plan"):
        assert section in doc
    assert "Reisman" in doc and "IRB PROTOCOL NUMBER" in doc


def test_protocol_expected_values_live():
    """The expected-values table must reflect the configured agent/doses."""
    doc = generate_protocol(ProtocolConfig(agent="metoclopramide",
                                           start_dose_mg_day=10,
                                           max_dose_mg_day=20))
    assert "| 10 |" in doc and "| 20 |" in doc


def test_protocol_domperidone_flags_not_approved():
    doc = generate_protocol(ProtocolConfig(agent="domperidone"))
    assert "NOT FDA-approved" in doc
    assert "QT" in doc


def test_protocol_deterministic():
    a = generate_protocol(ProtocolConfig())
    b = generate_protocol(ProtocolConfig())
    assert a == b


# --- trans-masculine protocol skeleton ------------------------------------------------

def test_protocol_t_contains_required_sections():
    doc = generate_androgen_protocol(AndrogenProtocolConfig(
        site="Test Clinic", pi="Dr. X", expected_values_n=20))
    for section in ("Background", "Specific aims", "Design", "Population",
                    "Interventions", "Expected values", "Safety and stopping",
                    "Ethics and consent", "Statistical plan"):
        assert section in doc
    # The domain-specific instruments: trough-timed T and the CBC schedule.
    assert "trough" in doc.lower()
    assert "CBC" in doc
    assert "IRB PROTOCOL NUMBER" in doc
    assert "protocol-t" in doc  # signature line names its generator


def test_protocol_t_expected_values_live():
    """The expected-values table must reflect the configured doses."""
    doc = generate_androgen_protocol(AndrogenProtocolConfig(
        start_dose_mg=50, max_dose_mg=150, expected_values_n=20))
    assert "| 50 mg / 7 d |" in doc and "| 150 mg / 7 d |" in doc


def test_protocol_t_expected_values_bands_finite():
    """Bands over points: the table must carry p10/p50/p90 triples, and the
    trough band must sit below the peak band (peak-trough structure visible)."""
    import re
    doc = generate_androgen_protocol(AndrogenProtocolConfig(expected_values_n=30))
    rows = re.findall(r"^\| (\d+) mg / 7 d \| (\d+) / (\d+) / (\d+) \| "
                      r"(\d+) \| (\d+) \| (\d+) / (\d+) / (\d+) \|", doc, flags=re.M)
    assert len(rows) == 2
    for r in rows:
        t_p10, t_p50, t_p90, mean50, peak50 = (int(r[i]) for i in range(1, 6))
        assert t_p10 <= t_p50 <= t_p90 < peak50
        assert mean50 < peak50


def test_protocol_t_gel_has_no_trough_language():
    doc = generate_androgen_protocol(AndrogenProtocolConfig(
        agent="testosterone_gel", expected_values_n=20))
    assert "no meaningful trough" in doc or "no trough exists" in doc


def test_protocol_t_deterministic():
    a = generate_androgen_protocol(AndrogenProtocolConfig(expected_values_n=20))
    b = generate_androgen_protocol(AndrogenProtocolConfig(expected_values_n=20))
    assert a == b
