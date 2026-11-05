"""Participant CSV ingest: pumping/lab logs -> personalized inference.

Turns the spreadsheet a participant actually keeps into the analysis layer.
Designed for the titration study's data collection — one row per visit day.

CSV schema (header required; column order free):

    date,drop_date,weeks_post_drop,dose_mg_day,prolactin_ng_ml,milk_ml_day,notes

* ``weeks_post_drop`` — preferred time axis; if absent, computed from
  ``date`` and ``drop_date`` (ISO dates, 7-day weeks).
* ``dose_mg_day`` — galactagogue dose in effect that day (the titration
  signal). Required unless a constant dose is supplied by the caller.
* ``prolactin_ng_ml`` / ``milk_ml_day`` — at least one per row; rows with
  neither are skipped (counted, not fatal — partial diaries are normal).
* ``milk_ml_day`` may be recorded in fl oz by community habit; pass
  ``milk_unit="oz"`` to convert (1 fl oz = 29.5735 mL).

Fail-closed on malformed data: every bad row is collected and reported with
its line number before anything runs — a study log with typos should produce
a fixable error list, not silent garbage.

The observations feed :func:`hormone_sim.personalize` unchanged, so the
posterior mechanics, ESS diagnostics, and caveats there apply here too.
"""

from __future__ import annotations

import csv
from datetime import date, datetime
from pathlib import Path
from typing import List, Optional, Tuple

from .personalization import Observation, PersonalizedResult, personalize, posterior_predictive

FL_OZ_ML = 29.5735


class IngestError(ValueError):
    """Raised when the CSV has problems; ``.problems`` lists row-level fixes."""


def _parse_float(text: Optional[str], field: str, rowno: int,
                 problems: List[str]) -> Optional[float]:
    if text is None or text.strip() == "":
        return None
    try:
        return float(text)
    except ValueError:
        problems.append(f"row {rowno}: {field}={text!r} is not a number")
        return None


def _parse_date(text: Optional[str], field: str, rowno: int,
                problems: List[str]) -> Optional[date]:
    if text is None or text.strip() == "":
        return None
    for fmt in ("%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"):
        try:
            return datetime.strptime(text.strip(), fmt).date()
        except ValueError:
            continue
    problems.append(f"row {rowno}: {field}={text!r} is not a date "
                    "(use YYYY-MM-DD)")
    return None


def load_observations(
    path: str | Path,
    milk_unit: str = "ml",
    default_dose_mg_day: Optional[float] = None,
) -> Tuple[List[Observation], int]:
    """Read a participant log. Returns (observations, skipped_row_count).

    Raises :class:`IngestError` listing every problem row (fail closed).
    """
    milk_unit = milk_unit.lower()
    if milk_unit not in ("ml", "oz"):
        raise ValueError(f"milk_unit must be 'ml' or 'oz', got {milk_unit!r}")

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    observations: List[Observation] = []
    skipped = 0
    problems: List[str] = []

    with p.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for rowno, row in enumerate(reader, start=2):
            get = lambda k: (row.get(k) or "").strip()  # noqa: E731
            if not any(get(k) for k in row):
                continue  # blank line

            # Time axis: weeks_post_drop, or date - drop_date.
            weeks = _parse_float(get("weeks_post_drop"), "weeks_post_drop",
                                 rowno, problems)
            if weeks is None:
                d = _parse_date(get("date"), "date", rowno, problems)
                dd = _parse_date(get("drop_date"), "drop_date", rowno, problems)
                if d is not None and dd is not None:
                    weeks = (d - dd).days / 7.0
            if weeks is None:
                problems.append(f"row {rowno}: no time axis (need "
                                "weeks_post_drop, or date + drop_date)")
                continue
            if weeks < 0:
                problems.append(f"row {rowno}: weeks_post_drop={weeks:g} < 0")
                continue

            dose = _parse_float(get("dose_mg_day"), "dose_mg_day", rowno, problems)
            if dose is None:
                dose = default_dose_mg_day
            if dose is None or dose <= 0:
                problems.append(f"row {rowno}: missing dose_mg_day and no "
                                "default dose supplied")
                continue

            prl = _parse_float(get("prolactin_ng_ml"), "prolactin_ng_ml",
                               rowno, problems)
            milk = _parse_float(get("milk_ml_day"), "milk_ml_day", rowno, problems)
            if milk is not None and milk_unit == "oz":
                milk *= FL_OZ_ML
            if milk is not None and milk < 0:
                problems.append(f"row {rowno}: negative milk volume")
                continue

            if prl is None and milk is None:
                skipped += 1
                continue

            observations.append(Observation(
                weeks_after_drop=weeks,
                dose_mg_day=dose,
                prolactin_ng_ml=prl,
                milk_ml_day=milk,
            ))

    if problems:
        raise IngestError(
            "CSV ingest found problems:\n  " + "\n  ".join(problems)
        )
    if not observations:
        raise IngestError("CSV ingest found no usable observations.")
    return observations, skipped


def personalize_from_csv(
    agent: str,
    path: str | Path,
    default_dose_mg_day: Optional[float] = None,
    milk_unit: str = "ml",
    prep_weeks: float = 30.0,
    pumping_sessions_per_day: float = 8.0,
    n: int = 400,
    seed: int = 7,
) -> PersonalizedResult:
    """Load a participant log and run the personal inference end-to-end."""
    obs, skipped = load_observations(path, milk_unit=milk_unit,
                                     default_dose_mg_day=default_dose_mg_day)
    return personalize(
        agent=agent, observations=obs, prep_weeks=prep_weeks,
        pumping_sessions_per_day=pumping_sessions_per_day,
        n=n, seed=seed,
    )


def plot_personal(
    result: PersonalizedResult,
    observations: List[Observation],
    prior_bands: Optional[Tuple] = None,
    out_path: Optional[str] = None,
    prep_weeks: float = 30.0,
    title: str = "Personalized supply posterior",
):
    """Plot the weighted posterior supply bands with observation overlays."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lo, med, hi = result.band(0.8)
    weeks = result.weeks - prep_weeks  # x-axis: weeks after the drop

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    if prior_bands is not None:
        plo, pmed, phi = prior_bands
        ax.fill_between(weeks, plo, phi, color="#999999", alpha=0.18,
                        label="population prior (80%)")
        ax.plot(weeks, pmed, color="#999999", lw=1.2, ls="--")
    ax.fill_between(weeks, lo, hi, color="#2980b9", alpha=0.30,
                    label="personal posterior (80%)")
    ax.plot(weeks, med, color="#2980b9", lw=2, label="posterior median")

    prl_obs = [o for o in observations if o.prolactin_ng_ml is not None]
    milk_obs = [o for o in observations if o.milk_ml_day is not None]
    if milk_obs:
        ax.scatter([o.weeks_after_drop for o in milk_obs],
                   [o.milk_ml_day for o in milk_obs],
                   color="#c0392b", zorder=5, label="measured milk")
    if prl_obs:
        # Prolactin observations are on a different scale; mark on the axis
        # top as reference ticks.
        for o in prl_obs:
            ax.axvline(o.weeks_after_drop, color="#27ae60", lw=0.8, alpha=0.5)
        ax.plot([], [], color="#27ae60", lw=0.8, alpha=0.5,
                label="prolactin lab draws")

    ax.set_xlabel("Weeks after the drop")
    ax.set_ylabel("Milk (mL/day)")
    ax.set_title(f"{title} — ESS {result.ess:.0f}/{result.weights.size}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=130)
        plt.close(fig)
    return fig
