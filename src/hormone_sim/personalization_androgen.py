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
from dataclasses import dataclass, field, replace
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
from .measurement import (
    BodySample,
    DomainPosterior,
    MeasurementDomain,
    DomainSweepRow,
    personalize_domain,
    posterior_predictive_table,
    population_sweep,
)


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


# params row layout shared by the domain's likelihood and predictive table.
_ANDROGEN_PARAM_NAMES = ["half_life", "vd", "ec50", "base_t", "base_e2",
                         "hct_rise", "hct_tau", "hct_baseline"]


@dataclass
class AndrogenDomain(MeasurementDomain):
    """The trans-masculine measurement model (see
    :mod:`hormone_sim.measurement`).

    Sampled body: T PK (half-life, Vd), HPG drive sensitivity (EC50),
    endogenous baselines, and the Hct response (rise ceiling, tau,
    baseline). Measurements: serum testosterone (index into the body's PK
    curve at the logged dose) and hematocrit (separate Hct-model
    simulation). The same descriptors also drive the posterior-predictive
    table and the population sweep — the T domain's titration language
    (trough/mean/peak, trough-in-range fraction) is presentation, the
    machinery is shared.
    """

    agent: str
    body_weight_kg: float = 70.0
    baseline_t_ngdl: Optional[float] = None   # intersex packs override here
    baseline_e2_pgml: Optional[float] = None
    dt_hours: float = 2.0
    total_weeks: Optional[float] = None       # personalization horizon
    pred_weeks: float = 12.0                  # predictive-table horizon
    trough_range_ngdl: Tuple[float, float] = (300.0, 1000.0)  # sweep criterion

    name: str = "androgen"
    param_names: List[str] = field(
        default_factory=lambda: list(_ANDROGEN_PARAM_NAMES))

    # predictive-table presentation (shared renderer in measurement.py)
    predictive_title: str = "Posterior predictive — T exposure at candidate regimens:"
    predictive_candidate_header: str = f"{'dose mg':>8} | {'q d':>4}"
    predictive_metrics: List[Tuple[str, int]] = field(default_factory=lambda: [
        ("trough p50", 10), ("mean p50", 8), ("peak p50", 8), ("Hct plat.", 8)])
    predictive_suffix: str = "(weighted medians over YOUR posterior)"
    predictive_rule: str = "-" * 78
    predictive_notes: Tuple[str, ...] = (
        "  Titration reading: weekly vs biweekly at the same weekly",
        "  dose buys the same mean with different peaks/troughs — pick",
        "  the swing that matches how your body handles troughs, and",
        "  check the Hct column against the monitoring thresholds.",
    )

    def __post_init__(self) -> None:
        self.a = get_androgen(self.agent)
        self._base_t = (BASELINE_T_NGDL if self.baseline_t_ngdl is None
                        else float(self.baseline_t_ngdl))
        self._base_e2 = (BASELINE_E2_PGML if self.baseline_e2_pgml is None
                         else float(self.baseline_e2_pgml))

    # --- personalization -----------------------------------------------------
    def sample_body(self, rng: np.random.Generator,
                    observations) -> BodySample:
        if self.total_weeks is None and observations:
            # Resolved lazily from the log when the caller didn't set it
            # (same rule as the public wrapper: span + 4 weeks).
            self.total_weeks = max(o.weeks_on_t for o in observations) + 4.0
        half_life = max(0.5, 2.0 * rng.lognormal(0.0, 0.30))
        vd = max(0.5, 2.0 * rng.lognormal(0.0, 0.30))
        ec50 = max(20.0, 170.0 * rng.lognormal(0.0, 0.30))
        b_t = max(1.0, self._base_t * rng.lognormal(0.0, 0.30))
        b_e2 = max(1.0, self._base_e2 * rng.lognormal(0.0, 0.30))
        hct_rise = max(1.0, 13.0 * rng.lognormal(0.0, 0.30))
        hct_tau = max(2.0, 10.0 * rng.lognormal(0.0, 0.40))
        hct_base = float(np.clip(rng.normal(40.0, 2.0), 33.0, 47.0))
        return BodySample(params=np.array(
            [half_life, vd, ec50, b_t, b_e2, hct_rise, hct_tau, hct_base]))

    def log_likelihood(self, body: BodySample, obs) -> float:
        if obs.testosterone_ngdl is None and obs.hct_pct is None:
            return 0.0
        half_life, vd, ec50, b_t, b_e2, hct_rise, hct_tau, hct_base = body.params
        lw = 0.0
        # The body's T curve over the full observation span, simulated at the
        # first measured observation's regimen (steady-state per-dose
        # approximation, as documented). Sampled Vd rides on a copy of the
        # agent (same depot kinetics, individual volume).
        if "t" not in body.curves:
            a_i = replace(self.a, vd_l_kg=vd)
            pk = simulate_androgen(
                a_i, obs.dose_mg, obs.interval_days, self.total_weeks,
                dt_hours=self.dt_hours, e2_ec50_ngdl=ec50,
                body_weight_kg=self.body_weight_kg,
                baseline_t_ngdl=b_t, baseline_e2_pgml=b_e2,
            )
            body.curves["t"] = pk["t_ngdl"]
            body.weeks = pk["hours"] / (7 * 24.0)
        if obs.testosterone_ngdl is not None:
            t_curve = body.curves["t"]
            idx = min(int(round(obs.weeks_on_t * 7 * 24 / self.dt_hours)),
                      t_curve.size - 1)
            pred = float(t_curve[idx])
            sigma = max(0.20 * obs.testosterone_ngdl, 30.0)
            lw += -0.5 * ((pred - obs.testosterone_ngdl) / sigma) ** 2
        if obs.hct_pct is not None:
            hm = HctModel(baseline_hct=hct_base, max_rise_points=hct_rise,
                          tau_wk=hct_tau)
            hr = simulate_hct(self.a, obs.dose_mg, obs.interval_days,
                              weeks=max(obs.weeks_on_t + 1.0, 2.0),
                              model=hm, body_weight_kg=self.body_weight_kg)
            idx2 = min(int(round(obs.weeks_on_t / 0.25)), hr.hct.size - 1)
            pred_h = float(hr.hct[idx2])
            sigma_h = max(0.05 * obs.hct_pct, 1.5)
            lw += -0.5 * ((pred_h - obs.hct_pct) / sigma_h) ** 2
        return lw

    def finalize(self, bodies: List[BodySample], observations) -> None:
        """Hct curves for every body at the participant's dominant regimen,
        aligned to the T curve's week axis."""
        dominant = max(observations, key=lambda o: o.weeks_on_t)
        weeks_axis = next((b.weeks for b in bodies if b.weeks is not None),
                          None)
        for b in bodies:
            half_life, vd, ec50, b_t, b_e2, hct_rise, hct_tau, hct_base = b.params
            hm = HctModel(baseline_hct=hct_base, max_rise_points=hct_rise,
                          tau_wk=hct_tau)
            hr = simulate_hct(self.a, dominant.dose_mg, dominant.interval_days,
                              weeks=self.total_weeks, model=hm,
                              body_weight_kg=self.body_weight_kg)
            if weeks_axis is None:
                # No body has a measured T observation (unreachable via the
                # CSV loaders, which require at least one measurement):
                # anchor the axes on the Hct simulation itself.
                weeks_axis = hr.weeks
            if "t" not in b.curves:
                b.curves["t"] = np.zeros(weeks_axis.size)
            b.curves["hct"] = np.interp(weeks_axis, hr.weeks, hr.hct)
            if b.weeks is None:
                b.weeks = weeks_axis

    def curve_names(self) -> Tuple[str, ...]:
        return ("t", "hct")

    def population_curve_names(self) -> Tuple[str, ...]:
        """Population sweeps ride the PK only (the Hct plateau column lives
        in the protocol generator's expected-values table)."""
        return ("t",)

    # --- posterior-predictive table ------------------------------------------
    def candidate_label(self, candidate: Tuple[float, float]) -> str:
        dose, interval = candidate
        return f"{dose:>8.1f} | {interval:>4.1f}"

    def predict_candidate(self, params: np.ndarray,
                          candidate: Tuple[float, float]) -> List[float]:
        dose, interval = candidate
        half_life, vd, ec50, b_t, b_e2, hct_rise, hct_tau, hct_base = params
        a_i = replace(self.a, vd_l_kg=vd)
        pk = simulate_androgen(
            a_i, dose, interval, self.pred_weeks, dt_hours=self.dt_hours,
            e2_ec50_ngdl=ec50, body_weight_kg=self.body_weight_kg,
            baseline_t_ngdl=b_t, baseline_e2_pgml=b_e2,
        )
        ss = pk["t_ngdl"][pk["t_ngdl"].size // 2:]
        hm = HctModel(baseline_hct=hct_base, max_rise_points=hct_rise,
                      tau_wk=hct_tau)
        hr = simulate_hct(self.a, dose, interval, weeks=self.pred_weeks,
                          model=hm, body_weight_kg=self.body_weight_kg)
        return [float(ss.min()), float(ss.mean()), float(ss.max()),
                hr.steady_state_hct()]

    # --- population sweep ------------------------------------------------------
    def population_body(self, rng: np.random.Generator,
                        candidate: Tuple[float, float]) -> BodySample:
        dose, interval = candidate
        half_life = max(0.5, self.a.half_life_h * rng.lognormal(0.0, 0.30))
        vd = max(0.5, self.a.vd_l_kg * rng.lognormal(0.0, 0.30))
        ec50 = max(20.0, 170.0 * rng.lognormal(0.0, 0.30))
        b_t = max(1.0, self._base_t * rng.lognormal(0.0, 0.30))
        b_e2 = max(1.0, self._base_e2 * rng.lognormal(0.0, 0.30))
        weeks_pk = max(8.0, 2.0 * interval / 7.0)
        a_i = replace(self.a, vd_l_kg=vd)
        pk = simulate_androgen(
            a_i, dose, interval, weeks_pk, dt_hours=self.dt_hours,
            e2_ec50_ngdl=ec50, body_weight_kg=self.body_weight_kg,
            baseline_t_ngdl=b_t, baseline_e2_pgml=b_e2,
        )
        body = BodySample(
            params=np.array([half_life, vd, ec50, b_t, b_e2,
                             13.0, 10.0, 40.0]),
            weeks=pk["hours"] / (7 * 24.0),
        )
        body.curves["t"] = pk["t_ngdl"]
        return body

    def sweep_summary(self, curves: dict) -> dict:
        t = curves["t"]
        ss = t[t.size // 2:]
        return {"trough": float(ss.min()), "mean": float(ss.mean()),
                "peak": float(ss.max())}

    def sweep_hit(self, summary: dict) -> bool:
        lo, hi = self.trough_range_ngdl
        return lo <= summary["trough"] <= hi


@dataclass
class AndrogenPersonalizedResult(DomainPosterior):
    """Posterior over modeled bodies given a trans-masculine participant's log.

    The T view of :class:`~hormone_sim.measurement.DomainPosterior`.
    """

    @property
    def t_ngdl(self) -> np.ndarray:
        return self.curves["t"]

    @property
    def hct(self) -> np.ndarray:
        return self.curves["hct"]

    def band(self, q: float = 0.8, which: str = "t") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Weighted (lower, median, upper) curves for central mass ``q``.

        ``which`` selects the outcome: ``"t"`` (testosterone, ng/dL) or
        ``"hct"`` (hematocrit, %). Weighted percentiles of the ensemble ARE
        the posterior; the unweighted percentiles (see
        :func:`plot_personal_androgen`) are the population prior.
        """
        if which not in ("t", "hct"):
            raise ValueError(f"which must be 't' or 'hct', got {which!r}")
        return super().band(q, which=which)

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
    """Reweight the androgen ensemble by a participant's T labs and CBCs.

    The sampling, likelihood weighting, and ESS machinery are the study
    layer's shared engine
    (:func:`hormone_sim.measurement.personalize_domain`) driven by
    :class:`AndrogenDomain`.
    """
    if not observations:
        raise ValueError("personalize_androgen() needs at least one observation")

    max_week = max(o.weeks_on_t for o in observations)
    total_weeks = max_week + 4.0

    domain = AndrogenDomain(
        agent=agent, body_weight_kg=body_weight_kg,
        baseline_t_ngdl=baseline_t_ngdl, baseline_e2_pgml=baseline_e2_pgml,
        dt_hours=dt_hours, total_weeks=total_weeks,
    )
    post = personalize_domain(domain, observations, n=n, seed=seed)
    return AndrogenPersonalizedResult(
        domain=domain, observations=post.observations, weights=post.weights,
        ess=post.ess, param_names=post.param_names, params=post.params,
        curves=post.curves, weeks=post.weeks,
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
    (dose_mg, interval_days) regimens — the trans-masculine titration table,
    rendered by the domain-generic table machinery."""
    domain = AndrogenDomain(agent=agent, body_weight_kg=body_weight_kg,
                            dt_hours=dt_hours, pred_weeks=weeks)
    return posterior_predictive_table(domain, result, candidates)


# --- Population sweep over (dose, interval) regimens -----------------------------
def androgen_population_sweep(
    agent: str,
    candidates: List[Tuple[float, float]],
    n: int = 100,
    seed: int = 42,
    body_weight_kg: float = 70.0,
    baseline_t_ngdl: Optional[float] = None,
    baseline_e2_pgml: Optional[float] = None,
    trough_range_ngdl: Tuple[float, float] = (300.0, 1000.0),
    dt_hours: float = 2.0,
) -> List[DomainSweepRow]:
    """Sweep (dose_mg, interval_days) regimens; quantify between-body
    heterogeneity in trough/mean/peak exposure and **what fraction of
    modeled bodies land their trough inside the target range** — the
    trans-masculine analog of the galactagogue low-dose sweep, on the same
    domain-generic machinery (:func:`hormone_sim.measurement.population_sweep`).

    The same seed is used for every candidate, so rows are paired — the same
    "bodies" are evaluated at each regimen.
    """
    domain = AndrogenDomain(
        agent=agent, body_weight_kg=body_weight_kg,
        baseline_t_ngdl=baseline_t_ngdl, baseline_e2_pgml=baseline_e2_pgml,
        trough_range_ngdl=trough_range_ngdl, dt_hours=dt_hours,
    )
    return population_sweep(domain, candidates, n, seed)


def androgen_sweep_table(rows: List[DomainSweepRow],
                         agent: str,
                         trough_range_ngdl: Tuple[float, float] = (300.0, 1000.0)
                         ) -> str:
    """Human-readable regimen sweep table (the `sweep-t` output)."""
    lo, hi = trough_range_ngdl
    lines = [
        f"T regimen sweep — {agent} (paired samples; n per regimen = "
        f"same bodies at each regimen):",
        f"{'regimen':>16} | {'trough p10/p50/p90':>20} | {'mean p50':>8} | "
        f"{'peak p50':>8} | {f'% trough in {lo:.0f}-{hi:.0f}':>18}",
        "-" * 84,
    ]
    for r in rows:
        tr = r.metrics["trough"]
        lines.append(
            f"{r.label:>16} | {tr[0]:>6.0f}/{tr[1]:>5.0f}/{tr[2]:>5.0f} | "
            f"{r.metric('mean', 'p50'):>8.0f} | {r.metric('peak', 'p50'):>8.0f} | "
            f"{r.frac_hit * 100:>17.1f}%"
        )
    lines += [
        "",
        "  Reading this: weekly vs biweekly at the same weekly dose buys the",
        "  same mean with different swings — the peak/trough columns are the",
        "  scheduling lesson, and the % column is the heterogeneity question",
        "  (titration-to-effect, not fixed dosing).",
        "",
        "  Caveats: parameter CVs are judgment calls, not measured population",
        "  spread; bands are population-typical estimation, never dosing advice;",
        "  no controlled dose-finding trials in trans-masculine populations.",
    ]
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


def plot_personal_androgen(
    result: AndrogenPersonalizedResult,
    observations: List[AndrogenObservation],
    out_path: Optional[str] = None,
    title: str = "Personalized T posterior",
):
    """Plot the weighted posterior T/Hct bands with observation overlays.

    The trans-masculine mirror of :func:`hormone_sim.ingest.plot_personal`.
    Top panel: serum testosterone with the male reference range (300-1000
    ng/dL, Endocrine Society CPG, trough-adjusted). Bottom panel: hematocrit
    with the commonly used evaluation/hold thresholds — markers, not
    decisions (see :mod:`hormone_sim.erythrocytosis` for the framing).

    The "population prior" band here is the SAME ensemble with the
    importance weights removed: every modeled body was simulated at the
    logged dose, so the unweighted percentiles are exactly the population
    prior at that regimen — a like-for-like posterior-vs-prior comparison,
    no extra simulation needed.

    Labs are assumed steady-state for the logged dose (the personalization
    caveat); the figure shows the posterior, the caveats live in
    :meth:`AndrogenPersonalizedResult.report`.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def unweighted_band(matrix: np.ndarray, q: float = 0.8):
        lo = np.percentile(matrix, (1.0 - q) / 2.0 * 100.0, axis=0)
        hi = np.percentile(matrix, (1.0 + q) / 2.0 * 100.0, axis=0)
        med = np.median(matrix, axis=0)
        return lo, med, hi

    weeks = result.weeks
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True,
                                   constrained_layout=True)

    # --- T panel -----------------------------------------------------------
    p_lo, p_med, p_hi = unweighted_band(result.t_ngdl)
    ax1.fill_between(weeks, p_lo, p_hi, color="#999999", alpha=0.18,
                     label="population prior (80%)")
    ax1.plot(weeks, p_med, color="#999999", lw=1.2, ls="--")
    t_lo, t_med, t_hi = result.band(0.8, which="t")
    ax1.fill_between(weeks, t_lo, t_hi, color="#2980b9", alpha=0.30,
                     label="personal posterior (80%)")
    ax1.plot(weeks, t_med, color="#2980b9", lw=2, label="posterior median")
    ax1.axhspan(300, 1000, color="#16a085", alpha=0.08,
                label="male reference 300-1000 (trough-adjusted)")
    t_obs = [o for o in observations if o.testosterone_ngdl is not None]
    if t_obs:
        ax1.scatter([o.weeks_on_t for o in t_obs],
                    [o.testosterone_ngdl for o in t_obs],
                    color="#c0392b", zorder=5, label="measured T")
    ax1.set_ylabel("Testosterone (ng/dL)")
    ax1.legend(fontsize=7, loc="upper right")
    ax1.grid(alpha=0.3)

    # --- Hct panel -----------------------------------------------------------
    hp_lo, hp_med, hp_hi = unweighted_band(result.hct)
    ax2.fill_between(weeks, hp_lo, hp_hi, color="#999999", alpha=0.18,
                     label="population prior (80%)")
    ax2.plot(weeks, hp_med, color="#999999", lw=1.2, ls="--")
    h_lo, h_med, h_hi = result.band(0.8, which="hct")
    ax2.fill_between(weeks, h_lo, h_hi, color="#2980b9", alpha=0.30,
                     label="personal posterior (80%)")
    ax2.plot(weeks, h_med, color="#2980b9", lw=2, label="posterior median")
    ax2.axhline(52, color="#000", ls=":", lw=0.8)
    ax2.axhline(54, color="#000", ls="--", lw=0.8)
    ax2.text(weeks[1], 52.15, "evaluate (practice)", fontsize=6, color="#000")
    ax2.text(weeks[1], 54.15, "hold/reduce zone (practice)", fontsize=6, color="#000")
    hct_obs = [o for o in observations if o.hct_pct is not None]
    if hct_obs:
        ax2.scatter([o.weeks_on_t for o in hct_obs],
                    [o.hct_pct for o in hct_obs],
                    color="#c0392b", zorder=5, label="measured Hct")
    ax2.set_ylabel("Hematocrit (%)")
    ax2.set_xlabel("Weeks on T")
    ax2.legend(fontsize=7, loc="lower right")
    ax2.grid(alpha=0.3)

    fig.suptitle(f"{title} — ESS {result.ess:.0f}/{result.weights.size}",
                 fontsize=11)
    if out_path:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=130)
    return fig
