"""Trans-masculine personalization: inference from T logs and CBCs.

The study-layer mirror of :mod:`hormone_sim.personalization` (which is
lactation-coupled). A trans-masculine participant's log records what's
actually measurable on T:

* dose (mg per injection) and interval (days) — or daily gel dose;
* serum total testosterone (ng/dL) — the direct titration signal;
* hematocrit (%) — the lab-monitorable safety measure.

Mechanism is the same importance-sampling scheme: sample each modeled body's
PK (half-life, volume), drive sensitivity (HPG EC50), endogenous baseline,
and Hct response parameters; weight by the likelihood of the observations;
report posterior medians, ESS, and a **posterior predictive over candidate
(dose, interval) regimens** — for T the interesting columns are *trough*,
*mean*, and *peak*, because schedule (weekly vs biweekly) changes the swing,
not just the level. That is the model's one genuinely useful scheduling
lesson, personalized.

Caveats carried over, plus two T-specific ones:

* Observations are assumed to be at steady state for the dose in effect
  (injection intervals are days-to-weeks; a lab drawn 3 days after a dose
  change is mid-transient — the model will under-weight it implicitly
  through the sigma, but the honest move is to log labs at steady state).
* Sampled per-body baselines include the intersex packs' territory (adrenal
  T in CAH, low testicular T in Klinefelter): pass the participant's pack
  via ``baseline_t_ngdl``/``baseline_e2_pgml`` when one applies.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from datetime import date, datetime
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from .androgens import (
    Androgen,
    get_androgen,
    simulate_androgen,
    BASELINE_T_NGDL,
    BASELINE_E2_PGML,
)
from .erythrocytosis import HctModel, simulate_hct
from .lactation import LactationProtocol  # noqa: F401  (type-compat only)


@dataclass
class AndrogenObservation:
    """One measured datapoint from a trans-masculine participant's log."""

    weeks_on_t: float
    dose_mg: float
    interval_days: float
    testosterone_ngdl: Optional[float] = None
    hct_pct: Optional[float] = None
    notes: str = ""

    def describe(self) -> str:
        bits = [f"wk {self.weeks_on_t:g} on T, "
                f"{self.dose_mg:g} mg/{self.interval_days:g} d"]
        if self.testosterone_ngdl is not None:
            bits.append(f"T {self.testosterone_ngdl:g} ng/dL")
        if self.hct_pct is not None:
            bits.append(f"Hct {self.hct_pct:g}%")
        return ", ".join(bits)


@dataclass
class AndrogenPersonalizedResult:
    """Posterior over modeled bodies given a trans-masculine participant's log."""

    observations: List[AndrogenObservation]
    weights: np.ndarray
    ess: float
    param_names: List[str]
    params: np.ndarray        # (n, k): half_life, vd, ec50, base_t, base_e2,
                              #          hct_rise, hct_tau, hct_baseline
    weeks: np.ndarray         # predictive time axis (weeks on T)
    t_ngdl: np.ndarray        # (n, t) unweighted T curves at the mean dose
    hct: np.ndarray           # (n, t) unweighted Hct curves

    def weighted_percentile(self, values: np.ndarray, p: float) -> float:
        idx = np.argsort(values)
        v = values[idx]
        w = self.weights[idx]
        cdf = np.cumsum(w)
        cdf /= cdf[-1]
        return float(np.interp(p / 100.0, cdf, v))

    def posterior_param(self, name: str, p: float = 50.0) -> float:
        return self.weighted_percentile(
            self.params[:, self.param_names.index(name)], p)

    def report(self) -> str:
        lines = [
            f"Personalized posterior (T) from {len(self.observations)} observation(s):",
        ]
        for o in self.observations:
            lines.append(f"  - {o.describe()}")
        lines += [
            f"  Effective sample size: {self.ess:.0f} of {self.weights.size}"
            + ("  (posterior ~= prior; the data hasn't moved the ensemble)"
               if self.ess > 0.8 * self.weights.size else
               "  (the observations are meaningfully constraining)"),
            "",
            "  Posterior medians (population prior -> your body):",
            f"    T half-life:        {self.posterior_param('half_life'):.1f} h",
            f"    T Vd:               {self.posterior_param('vd'):.1f} L/kg",
            f"    HPG drive EC50:     {self.posterior_param('ec50'):.0f} ng/dL",
            f"    endogenous T:       {self.posterior_param('base_t'):.0f} ng/dL",
            f"    Hct rise ceiling:   +{self.posterior_param('hct_rise'):.0f} points",
            "",
            "  Caveats: labs are assumed steady-state for the logged dose;",
            "  the posterior reweights a population prior (postpartum-cis-",
            "  adjacent literature for Hct; extrapolated throughout) — it does",
            "  not launder the prior's weaknesses. Use under clinical care.",
        ]
        return "\n".join(lines)


def personalize_androgen(
    agent: str,
    observations: List[AndrogenObservation],
    n: int = 200,
    seed: int = 7,
    body_weight_kg: float = 70.0,
    baseline_t_ngdl: Optional[float] = None,
    baseline_e2_pgml: Optional[float] = None,
    dt_hours: float = 2.0,
) -> AndrogenPersonalizedResult:
    """Reweight the androgen ensemble by a participant's T labs and CBCs."""
    if not observations:
        raise ValueError("personalize_androgen() needs at least one observation")
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
    hct_curves = None
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
        # One T curve per body over the full observation span (steady-state
        # per-dose approximation, as documented). The sampled Vd rides on a
        # copy of the agent (same depot kinetics, individual volume).
        t_curve = None
        weeks_i = None
        a_i = replace(a, vd_l_kg=vd)
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
            # T likelihood.
            if o.testosterone_ngdl is not None:
                idx = min(int(round(o.weeks_on_t * 7 * 24 / dt_hours)),
                          t_curve.size - 1)
                pred = float(t_curve[idx])
                sigma = max(0.20 * o.testosterone_ngdl, 30.0)
                lw += -0.5 * ((pred - o.testosterone_ngdl) / sigma) ** 2
            # Hct likelihood (separate sampled model).
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
            hct_curves = np.empty((n, t_curve.size))
            weeks_axis = weeks_i
        if t_curve is not None:
            t_curves[i] = t_curve
        else:
            t_curves[i] = 0.0

    # Hct curves: computed per body at the participant's dominant regimen,
    # aligned to the T curve's week axis.
    dominant = max(observations, key=lambda o: o.weeks_on_t)
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

    return AndrogenPersonalizedResult(
        observations=list(observations), weights=w, ess=float(ess),
        param_names=param_names, params=P, weeks=weeks_axis,
        t_ngdl=t_curves, hct=hct_curves,
    )


def posterior_predictive_androgen(
    result: AndrogenPersonalizedResult,
    agent: str,
    candidates: List[Tuple[float, float]],
    weeks: float = 12.0,
    body_weight_kg: float = 70.0,
    dt_hours: float = 2.0,
) -> str:
    """Weighted trough/mean/peak + Hct-plateau predictions at candidate
    (dose_mg, interval_days) regimens — the trans-masculine titration table."""
    a = get_androgen(agent)
    names = result.param_names
    lines = [
        "Posterior predictive — T exposure at candidate regimens:",
        f"{'dose mg':>8} | {'q d':>4} | {'trough p50':>10} | {'mean p50':>8} | "
        f"{'peak p50':>8} | {'Hct plat.':>8}   (weighted medians over YOUR posterior)",
        "-" * 78,
    ]
    for dose, interval in candidates:
        troughs, means, peaks, hcts = [], [], [], []
        for i in range(result.weights.size):
            a_i = replace(a, vd_l_kg=result.params[i, names.index("vd")])
            pk = simulate_androgen(
                a_i, dose, interval, weeks, dt_hours=dt_hours,
                e2_ec50_ngdl=result.params[i, names.index("ec50")],
                body_weight_kg=body_weight_kg,
                baseline_t_ngdl=result.params[i, names.index("base_t")],
                baseline_e2_pgml=result.params[i, names.index("base_e2")],
            )
            ss = pk["t_ngdl"][pk["t_ngdl"].size // 2:]
            troughs.append(ss.min())
            means.append(ss.mean())
            peaks.append(ss.max())
            hm = HctModel(baseline_hct=result.params[i, names.index("hct_baseline")],
                          max_rise_points=result.params[i, names.index("hct_rise")],
                          tau_wk=result.params[i, names.index("hct_tau")])
            hr = simulate_hct(a, dose, interval, weeks=weeks, model=hm,
                              body_weight_kg=body_weight_kg)
            hcts.append(hr.steady_state_hct())
        troughs, means, peaks, hcts = map(np.array, (troughs, means, peaks, hcts))
        lines.append(
            f"{dose:>8.1f} | {interval:>4.1f} | "
            f"{result.weighted_percentile(troughs, 50):>10.0f} | "
            f"{result.weighted_percentile(means, 50):>8.0f} | "
            f"{result.weighted_percentile(peaks, 50):>8.0f} | "
            f"{result.weighted_percentile(hcts, 50):>8.0f}"
        )
    lines.append("")
    lines.append("  Titration reading: weekly vs biweekly at the same weekly")
    lines.append("  dose buys the same mean with different peaks/troughs — pick")
    lines.append("  the swing that matches how your body handles troughs, and")
    lines.append("  check the Hct column against the monitoring thresholds.")
    return "\n".join(lines)


# --- CSV ingest for T logs ------------------------------------------------------
def load_androgen_observations(
    path: str | Path,
    default_interval_days: Optional[float] = None,
) -> Tuple[List[AndrogenObservation], int]:
    """Read a trans-masculine participant log (fail-closed, like the trans-fem one).

    Schema (header required; order free):

        date,start_date,weeks_on_t,dose_mg,interval_days,testosterone_ngdl,hct_pct,notes

    ``weeks_on_t`` preferred; else computed from ``date``/``start_date``.
    Rows need at least one of testosterone_ngdl / hct_pct.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    observations: List[AndrogenObservation] = []
    skipped = 0
    problems: List[str] = []

    def _f(text, field, rowno):
        if text is None or str(text).strip() == "":
            return None
        try:
            return float(text)
        except ValueError:
            problems.append(f"row {rowno}: {field}={text!r} is not a number")
            return None

    def _d(text, field, rowno):
        if text is None or str(text).strip() == "":
            return None
        for fmt in ("%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"):
            try:
                return datetime.strptime(str(text).strip(), fmt).date()
            except ValueError:
                continue
        problems.append(f"row {rowno}: {field}={text!r} is not a date")
        return None

    with p.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for rowno, row in enumerate(reader, start=2):
            get = lambda k: (row.get(k) or "").strip()  # noqa: E731
            if not any(get(k) for k in row):
                continue
            weeks = _f(get("weeks_on_t"), "weeks_on_t", rowno)
            if weeks is None:
                d = _d(get("date"), "date", rowno)
                sd = _d(get("start_date"), "start_date", rowno)
                if d is not None and sd is not None:
                    weeks = (d - sd).days / 7.0
            if weeks is None:
                problems.append(f"row {rowno}: no time axis")
                continue
            if weeks < 0:
                problems.append(f"row {rowno}: weeks_on_t < 0")
                continue
            dose = _f(get("dose_mg"), "dose_mg", rowno)
            if dose is None or dose <= 0:
                problems.append(f"row {rowno}: missing/invalid dose_mg")
                continue
            interval = _f(get("interval_days"), "interval_days", rowno)
            if interval is None:
                interval = default_interval_days
            if interval is None or interval <= 0:
                problems.append(f"row {rowno}: missing interval_days and no default")
                continue
            t = _f(get("testosterone_ngdl"), "testosterone_ngdl", rowno)
            hct = _f(get("hct_pct"), "hct_pct", rowno)
            if t is None and hct is None:
                skipped += 1
                continue
            observations.append(AndrogenObservation(
                weeks_on_t=weeks, dose_mg=dose, interval_days=interval,
                testosterone_ngdl=t, hct_pct=hct,
                notes=get("notes"),
            ))

    if problems:
        from .ingest import IngestError
        raise IngestError("T-log ingest found problems:\n  "
                          + "\n  ".join(problems))
    if not observations:
        raise IngestError("T-log ingest found no usable observations.")
    return observations, skipped
