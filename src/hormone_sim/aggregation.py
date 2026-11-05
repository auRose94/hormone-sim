"""Cohort aggregation: de-identified multi-participant case-series summaries.

The running-case-series view the protocol's statistical plan (§9) consumes:
many participants' logs, one de-identified table. This is the study layer's
aggregation step — it does not re-run the model and it performs no
inference; it summarizes what was *observed*, per participant and across the
cohort, so a clinician can see a case series accumulate.

**Privacy shape (non-negotiable):**

* Participants are keyed by ``study_id`` — a code the site assigns. The
  manifest never carries names, and the report says so; study codes must
  not be names.
* The underlying logs remain identifiable health data (they carry dates and
  dosing histories). This module only produces the de-identified aggregate
  view; storage policy is the site's (protocol §10, AGENTS.md rule 7).

**Input: a cohort manifest CSV** (header required; order free):

    study_id,csv,domain

* ``study_id`` — de-identified code (unique; duplicates fail closed).
* ``csv`` — path to that participant's log, in the same schema the
  ``ingest`` (lactation) or ``ingest-t`` (T log) commands read. Paths are
  used as given; if not found, they are resolved relative to the manifest's
  own directory so a clinic can keep manifest + logs in one folder.
* ``domain`` — ``lactation`` or ``t_log``.

Fail-closed, like every ingest path in this package: problems are collected
per participant (with the study_id attached) and raised together before
anything is summarized.

One flat summary row shape serves both domains — measures a domain doesn't
produce are simply None in its rows. The uniform shape is deliberate: it is
the same table the statistical plan reads regardless of which domain a
participant belongs to.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from .ingest import IngestError, load_observations
from .personalization_androgen import load_androgen_observations

DOMAINS = ("lactation", "t_log")

SUMMARY_FIELDS = [
    "study_id", "domain", "n_observations", "week_start", "week_end",
    "dose_min", "dose_max",
    # lactation-specific (None for t_log rows)
    "n_prl", "prl_first", "prl_max",
    "n_milk", "milk_first", "milk_last", "weeks_to_milk_threshold",
    # t_log-specific (None for lactation rows)
    "interval_days", "n_t", "t_first", "t_last", "t_max",
    "n_hct", "hct_first", "hct_max",
]


@dataclass
class ParticipantSummary:
    """One participant's observed-data summary (no model output)."""

    study_id: str
    domain: str
    n_observations: int
    week_start: float
    week_end: float
    dose_min: float
    dose_max: float
    n_prl: Optional[int] = None
    prl_first: Optional[float] = None
    prl_max: Optional[float] = None
    n_milk: Optional[int] = None
    milk_first: Optional[float] = None
    milk_last: Optional[float] = None
    weeks_to_milk_threshold: Optional[float] = None
    interval_days: Optional[float] = None
    n_t: Optional[int] = None
    t_first: Optional[float] = None
    t_last: Optional[float] = None
    t_max: Optional[float] = None
    n_hct: Optional[int] = None
    hct_first: Optional[float] = None
    hct_max: Optional[float] = None

    def row(self) -> dict:
        """Flat dict in SUMMARY_FIELDS order (for the per-participant CSV)."""
        return {k: getattr(self, k) for k in SUMMARY_FIELDS}

    def describe(self) -> str:
        bits = [
            f"{self.study_id}: {self.n_observations} obs, "
            f"wk {self.week_start:g}-{self.week_end:g}, "
            f"dose {self.dose_min:g}-{self.dose_max:g}"
        ]
        if self.domain == "lactation":
            if self.prl_max is not None:
                bits.append(f"PRL {self.prl_first:g}->{self.prl_max:g} ng/mL")
            if self.milk_last is not None:
                bits.append(f"milk {self.milk_first:g}->{self.milk_last:g} mL/d")
            bits.append("50 mL/d at wk "
                        + (f"{self.weeks_to_milk_threshold:g}"
                           if self.weeks_to_milk_threshold is not None else "n/r"))
        else:
            if self.t_max is not None:
                bits.append(f"T {self.t_first:g}->{self.t_last:g} ng/dL "
                            f"(max {self.t_max:g})")
            if self.hct_max is not None:
                bits.append(f"Hct {self.hct_first:g}->{self.hct_max:g}%")
        return ", ".join(bits)


def summarize_lactation_log(study_id: str, obs,
                            milk_threshold_ml_day: float = 50.0
                            ) -> ParticipantSummary:
    """Summarize one participant's lactation log (observed data only)."""
    prl = [o.prolactin_ng_ml for o in obs if o.prolactin_ng_ml is not None]
    milk = [(o.weeks_after_drop, o.milk_ml_day) for o in obs
            if o.milk_ml_day is not None]
    weeks_to_thr = None
    if milk:
        reached = [w for w, m in milk if m >= milk_threshold_ml_day]
        if reached:
            weeks_to_thr = min(reached)
    return ParticipantSummary(
        study_id=study_id, domain="lactation",
        n_observations=len(obs),
        week_start=min(o.weeks_after_drop for o in obs),
        week_end=max(o.weeks_after_drop for o in obs),
        dose_min=min(o.dose_mg_day for o in obs),
        dose_max=max(o.dose_mg_day for o in obs),
        n_prl=len(prl),
        prl_first=(prl[0] if prl else None),
        prl_max=(max(prl) if prl else None),
        n_milk=len(milk),
        milk_first=(milk[0][1] if milk else None),
        milk_last=(milk[-1][1] if milk else None),
        weeks_to_milk_threshold=weeks_to_thr,
    )


def summarize_t_log(study_id: str, obs) -> ParticipantSummary:
    """Summarize one participant's T log (observed data only)."""
    ts = [(o.weeks_on_t, o.testosterone_ngdl) for o in obs
          if o.testosterone_ngdl is not None]
    hcts = [(o.weeks_on_t, o.hct_pct) for o in obs if o.hct_pct is not None]
    intervals = Counter(o.interval_days for o in obs)
    return ParticipantSummary(
        study_id=study_id, domain="t_log",
        n_observations=len(obs),
        week_start=min(o.weeks_on_t for o in obs),
        week_end=max(o.weeks_on_t for o in obs),
        dose_min=min(o.dose_mg for o in obs),
        dose_max=max(o.dose_mg for o in obs),
        interval_days=(intervals.most_common(1)[0][0] if intervals else None),
        n_t=len(ts),
        t_first=(ts[0][1] if ts else None),
        t_last=(ts[-1][1] if ts else None),
        t_max=(max(t for _, t in ts) if ts else None),
        n_hct=len(hcts),
        hct_first=(hcts[0][1] if hcts else None),
        hct_max=(max(h for _, h in hcts) if hcts else None),
    )


def _resolve_csv(path_text: str, manifest_dir: Path) -> Path:
    """As-given first (CWD-relative or absolute), then manifest-relative."""
    p = Path(path_text)
    if p.exists() or p.is_absolute():
        return p
    return manifest_dir / p


def load_cohort(
    manifest_path: str | Path,
    milk_unit: str = "ml",
    default_dose_mg_day: Optional[float] = None,
    default_interval_days: Optional[float] = None,
    milk_threshold_ml_day: float = 50.0,
) -> List[ParticipantSummary]:
    """Read a cohort manifest and summarize every participant, fail-closed."""
    p = Path(manifest_path)
    if not p.exists():
        raise FileNotFoundError(p)
    manifest_dir = p.parent

    summaries: List[ParticipantSummary] = []
    seen_ids: set = set()
    problems: List[str] = []

    with p.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for rowno, row in enumerate(reader, start=2):
            get = lambda k: (row.get(k) or "").strip()  # noqa: E731
            if not any(get(k) for k in row):
                continue
            study_id = get("study_id")
            csv_path = get("csv")
            domain = get("domain").lower()
            if not study_id or not csv_path or not domain:
                problems.append(f"manifest row {rowno}: need study_id, csv, "
                                "and domain")
                continue
            if study_id in seen_ids:
                problems.append(f"manifest row {rowno}: duplicate study_id "
                                f"{study_id!r}")
                continue
            seen_ids.add(study_id)
            if domain not in DOMAINS:
                problems.append(f"manifest row {rowno} ({study_id}): unknown "
                                f"domain {domain!r} (want {DOMAINS})")
                continue
            try:
                if domain == "lactation":
                    obs, _skipped = load_observations(
                        _resolve_csv(csv_path, manifest_dir),
                        milk_unit=milk_unit,
                        default_dose_mg_day=default_dose_mg_day)
                    summaries.append(summarize_lactation_log(
                        study_id, obs,
                        milk_threshold_ml_day=milk_threshold_ml_day))
                else:
                    obs, _skipped = load_androgen_observations(
                        _resolve_csv(csv_path, manifest_dir),
                        default_interval_days=default_interval_days)
                    summaries.append(summarize_t_log(study_id, obs))
            except (IngestError, FileNotFoundError) as e:
                problems.append(f"{study_id}: {e}")

    if problems:
        raise IngestError("Cohort ingest found problems:\n  "
                          + "\n  ".join(problems))
    if not summaries:
        raise IngestError("Cohort manifest produced no participant summaries.")
    return summaries


def cohort_rows(summaries: List[ParticipantSummary]) -> List[dict]:
    """Per-participant rows in SUMMARY_FIELDS order (the statistical plan's
    de-identified input table)."""
    return [s.row() for s in summaries]


def write_cohort_csv(summaries: List[ParticipantSummary], path: str | Path
                     ) -> None:
    """Write the per-participant de-identified table."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        w.writeheader()
        w.writerows(cohort_rows(summaries))


def _iqr(vals: List[float]) -> Optional[Tuple]:
    if len(vals) < 2:
        return None
    return (float(np.percentile(vals, 25)), float(np.percentile(vals, 75)))


def cohort_report(summaries: List[ParticipantSummary],
                  milk_threshold_ml_day: float = 50.0) -> str:
    """Human-readable de-identified cohort summary (the case-series view)."""
    lines = [
        f"Cohort aggregation — running case series "
        f"({len(summaries)} participant(s), "
        f"{sum(s.n_observations for s in summaries)} observation(s))",
        "",
    ]
    for domain in DOMAINS:
        part = [s for s in summaries if s.domain == domain]
        if not part:
            continue
        lines.append(f"  [{domain}]")
        for s in part:
            lines.append(f"    {s.describe()}")
        # Cohort-level aggregates over the per-participant summaries.
        if domain == "lactation":
            lasts = [s.milk_last for s in part if s.milk_last is not None]
            thr = [s for s in part
                   if s.weeks_to_milk_threshold is not None]
            iqr = _iqr(lasts)
            lines.append(
                f"    cohort (n={len(part)}): latest milk median "
                f"{np.median(lasts):.0f} mL/d"
                + (f" (IQR {iqr[0]:.0f}-{iqr[1]:.0f})" if iqr else "")
                + f"; {len(thr)}/{len(part)} reached "
                f"{milk_threshold_ml_day:g} mL/d"
                + (f" (median wk {np.median([s.weeks_to_milk_threshold for s in thr]):.1f})"
                   if thr else ""))
        else:
            lasts = [s.t_last for s in part if s.t_last is not None]
            hmax = [s.hct_max for s in part if s.hct_max is not None]
            hi = [s for s in part
                  if s.hct_max is not None and s.hct_max >= 52.0]
            iqr = _iqr(lasts)
            lines.append(
                f"    cohort (n={len(part)}): latest T median "
                f"{np.median(lasts):.0f} ng/dL"
                + (f" (IQR {iqr[0]:.0f}-{iqr[1]:.0f})" if iqr else "")
                + (f"; max Hct median {np.median(hmax):.1f}%" if hmax else "")
                + f"; {len(hi)}/{len(part)} with an Hct reading >= 52% "
                "(the evaluate threshold)")
        lines.append("")

    lines += [
        "  Privacy: keyed by study_id only; study codes must not be names.",
        "  Logs themselves remain identifiable health data — store per site",
        "  policy and report de-identified aggregates publicly.",
        "  Descriptive observed-data summary for the protocol's statistical",
        "  plan; no model inference here. Population model output remains",
        "  population-typical estimation, never dosing advice.",
    ]
    return "\n".join(lines)
