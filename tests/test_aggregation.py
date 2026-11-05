"""Tests for the cohort aggregation view (AGENTS.md §5a item 4).

The aggregation layer must be de-identified (study_id keyed only),
fail-closed (problems named per participant), and produce one uniform
summary-row shape across domains — the table the protocol's statistical
plan consumes.
"""

import csv

import pytest

from hormone_sim import (
    IngestError,
    load_cohort,
    cohort_report,
    cohort_rows,
    write_cohort_csv,
)
from hormone_sim.aggregation import SUMMARY_FIELDS


@pytest.fixture
def cohort_dir(tmp_path):
    """Manifest + one log per domain, paths relative to the manifest."""
    (tmp_path / "tf1.csv").write_text(
        "weeks_post_drop,dose_mg_day,prolactin_ng_ml,milk_ml_day\n"
        "0,10,58,\n"
        "2,10,72,95\n"
        "6,15,78,140\n")
    (tmp_path / "tm1.csv").write_text(
        "weeks_on_t,dose_mg,interval_days,testosterone_ngdl,hct_pct\n"
        "0,100,7,410,41.2\n"
        "6,100,7,880,44.0\n"
        "12,100,7,920,45.8\n")
    (tmp_path / "manifest.csv").write_text(
        "study_id,csv,domain\n"
        "TF-001,tf1.csv,lactation\n"
        "TM-001,tm1.csv,t_log\n")
    return tmp_path / "manifest.csv"


def test_cohort_loads_both_domains(cohort_dir):
    summaries = load_cohort(cohort_dir)
    assert [s.study_id for s in summaries] == ["TF-001", "TM-001"]
    tf, tm = summaries
    assert tf.domain == "lactation"
    assert tf.n_observations == 3 and tf.prl_max == 78.0
    assert tf.milk_last == 140.0
    assert tf.weeks_to_milk_threshold == pytest.approx(2.0)
    assert tm.domain == "t_log"
    assert tm.t_max == 920.0 and tm.hct_first == 41.2 and tm.hct_max == 45.8
    assert tm.interval_days == 7.0
    # Uniform row shape across domains: same fields, other domain's = None.
    assert set(cohort_rows(summaries)[0].keys()) == set(SUMMARY_FIELDS)
    assert tm.prl_max is None and tf.t_max is None


def test_cohort_report_content(cohort_dir):
    text = cohort_report(load_cohort(cohort_dir))
    assert "TF-001" in text and "TM-001" in text
    assert "study codes must not be names" in text
    assert "no model inference" in text
    assert "reached 50 mL/d" in text
    assert ">= 52%" in text  # the erythrocytosis evaluate threshold


def test_cohort_manifest_relative_paths(cohort_dir):
    """Bare filenames must resolve against the manifest's directory."""
    summaries = load_cohort(cohort_dir)
    assert len(summaries) == 2


def test_cohort_fails_closed_with_study_id(tmp_path):
    (tmp_path / "bad.csv").write_text(
        "weeks_post_drop,dose_mg_day,milk_ml_day\n"
        "2,,100\n")  # missing dose, no default -> problem
    (tmp_path / "manifest.csv").write_text(
        "study_id,csv,domain\nTF-009,bad.csv,lactation\n")
    with pytest.raises(IngestError) as e:
        load_cohort(tmp_path / "manifest.csv")
    assert "TF-009" in str(e.value)


def test_cohort_rejects_duplicate_and_unknown_domain(tmp_path):
    (tmp_path / "a.csv").write_text(
        "weeks_post_drop,dose_mg_day,milk_ml_day\n2,10,100\n")
    (tmp_path / "manifest.csv").write_text(
        "study_id,csv,domain\n"
        "X-1,a.csv,lactation\n"
        "X-1,a.csv,lactation\n"
        "X-2,a.csv,hrt\n")
    with pytest.raises(IngestError) as e:
        load_cohort(tmp_path / "manifest.csv")
    msg = str(e.value)
    assert "duplicate study_id" in msg and "unknown domain" in msg


def test_cohort_oz_passthrough(tmp_path):
    (tmp_path / "oz.csv").write_text(
        "weeks_post_drop,dose_mg_day,milk_ml_day\n2,10,3.3814\n")
    (tmp_path / "manifest.csv").write_text(
        "study_id,csv,domain\nTF-OZ,oz.csv,lactation\n")
    s = load_cohort(tmp_path / "manifest.csv", milk_unit="oz")[0]
    assert s.milk_last == pytest.approx(100.0, rel=0.01)


def test_cohort_missing_manifest():
    with pytest.raises(FileNotFoundError):
        load_cohort("/nonexistent/manifest.csv")


def test_cohort_csv_round_trip(cohort_dir, tmp_path):
    summaries = load_cohort(cohort_dir)
    out = tmp_path / "cohort.csv"
    write_cohort_csv(summaries, out)
    with out.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    assert rows[0]["study_id"] == "TF-001"
    assert rows[1]["t_max"] == "920.0"
    assert rows[0]["t_max"] == ""  # other domain's measure is blank, not 0


def test_participant_summary_no_model_output():
    """Guard the layer's honesty: the summary dataclass carries no model
    fields (no predictions, no posteriors) — observed data only."""
    model_words = ("weights", "ess", "posterior", "predicted", "band")
    for f in SUMMARY_FIELDS:
        assert not any(w in f for w in model_words)
