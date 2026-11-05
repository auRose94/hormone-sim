"""Participant-data ingest: personal inference from measured points.

This is the module that turns the model from a *population prior* into a
*personal inference tool* — the piece a titration study actually needs.

Mechanism (importance sampling, no MCMC):

1. Sample N modeled bodies: galactagogue Emax/ED50 + baseline prolactin
   (drug side) and the full lactation parameter set (gland side), each
   lognormal around population defaults.
2. For each measured observation — a prolactin level and/or a milk volume at
   a known week and dose — compute each body's prediction and weight it by
   the likelihood of the observation (normal, with pragmatic sigmas:
   ~15% + 3 ng/mL for prolactin assay + within-person variation; ~30% + 10
   mL/day for milk day-to-day variation).
3. The reweighted ensemble IS the posterior: weighted percentile bands of
   milk supply over time, posterior medians of the person's parameters,
   and a **posterior predictive for candidate next doses** — the titration
   decision ("if I go 10 → 15 mg/d, what supply should I expect?").

Diagnostics, honestly stated:

* **Effective sample size (ESS)** is reported. A few observations barely
  move the ensemble (posterior ≈ prior — the tool says so); many
  contradictory observations collapse the weights (ESS → small; the
  posterior then reflects the model's inability to fit, not certainty).
* Milk predictions assume the observation's dose was constant from the drop
  onward. Under real titration this overstates precision somewhat — treat
  the milk posterior as indicative; the *prolactin* posterior (steady-state,
  dose-linked) is the tighter, more trustworthy update.
* The posterior is only as good as the population prior it reweights — all
  the "no controlled trials in trans women" caveats from the rest of the
  package apply underneath it. What personalization adds is the *person's*
  data on top; it does not launder the prior's weaknesses.

Not medical advice: this is the analysis layer for a clinician-supervised
titration — the observations come from labs and pumping logs collected under
care.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import List, Optional, Tuple

import numpy as np

from .galactagogues import (
    get_galactagogue,
    PROLACTIN_BASELINE_NG_ML,
)
from .lactation import LactationProtocol, LactationModel, simulate_lactation
from .sensitivity import SensitivitySpec


@dataclass
class Observation:
    """One measured datapoint from a participant.

    ``weeks_after_drop`` is time since the hormone drop (the galactagogue
    start). ``dose_mg_day`` is the galactagogue dose in effect at that time.
    Either measurement may be None (e.g. a prolactin-only lab visit).
    """

    weeks_after_drop: float
    dose_mg_day: float
    prolactin_ng_ml: Optional[float] = None
    milk_ml_day: Optional[float] = None

    def describe(self) -> str:
        bits = [f"wk {self.weeks_after_drop:g} post-drop, {self.dose_mg_day:g} mg/d"]
        if self.prolactin_ng_ml is not None:
            bits.append(f"PRL {self.prolactin_ng_ml:g} ng/mL")
        if self.milk_ml_day is not None:
            bits.append(f"milk {self.milk_ml_day:g} mL/day")
        return ", ".join(bits)


@dataclass
class PersonalizedResult:
    """Posterior over modeled bodies given a participant's observations."""

    observations: List[Observation]
    weights: np.ndarray            # normalized importance weights, shape (n,)
    ess: float                     # effective sample size
    param_names: List[str]         # names in `params` columns
    params: np.ndarray             # shape (n, k) sampled parameter values
    weeks: np.ndarray              # time axis of predictive curves
    milk: np.ndarray               # shape (n, t) unweighted supply curves
    prior_spec: SensitivitySpec

    # --- weighted statistics -------------------------------------------------
    def weighted_percentile(self, values: np.ndarray, p: float) -> float:
        idx = np.argsort(values)
        v = values[idx]
        w = self.weights[idx]
        cdf = np.cumsum(w)
        cdf /= cdf[-1]
        return float(np.interp(p / 100.0, cdf, v))

    def band(self, q: float = 0.8) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Weighted (lower, median, upper) milk curves for central mass q."""
        lo_pct = (1.0 - q) / 2.0
        hi_pct = (1.0 + q) / 2.0
        t = self.milk.shape[1]
        lo = np.array([self.weighted_percentile(self.milk[:, j], lo_pct * 100)
                       for j in range(t)])
        med = np.array([self.weighted_percentile(self.milk[:, j], 50.0)
                        for j in range(t)])
        hi = np.array([self.weighted_percentile(self.milk[:, j], hi_pct * 100)
                       for j in range(t)])
        return lo, med, hi

    def posterior_param(self, name: str, p: float = 50.0) -> float:
        return self.weighted_percentile(self.params[:, self.param_names.index(name)], p)

    def report(self) -> str:
        lines = [
            f"Personalized posterior from {len(self.observations)} observation(s):",
        ]
        for o in self.observations:
            lines.append(f"  - {o.describe()}")
        lines += [
            f"  Effective sample size: {self.ess:.0f} of {self.weights.size}"
            + ("  (posterior ~= prior; the data has not moved the ensemble much)"
               if self.ess > 0.8 * self.weights.size else
               "  (the observations are meaningfully constraining)"),
            "",
            "  Posterior medians (population prior -> your body):",
            f"    prolactin Emax:      {self.posterior_param('emax'):.0f} ng/mL"
            f"  (prior median = population value)",
            f"    prolactin ED50:      {self.posterior_param('ed50'):.1f} mg/day",
            f"    baseline prolactin:  {self.posterior_param('baseline'):.1f} ng/mL",
            f"    milk Vmax:           {self.posterior_param('milk_vmax'):.0f} mL/day",
            f"    prolactin Km (milk): {self.posterior_param('milk_km'):.0f} ng/mL",
            "",
            "  Weighted supply bands at the CURRENT dose (mL/day, last 8 wk):",
            f"    p10 {self.weighted_percentile(self.milk[:, -1], 10):.0f} / "
            f"median {self.weighted_percentile(self.milk[:, -1], 50):.0f} / "
            f"p90 {self.weighted_percentile(self.milk[:, -1], 90):.0f}",
            "",
            "  Caveats: milk predictions assume each observation's dose was",
            "  constant from the drop (titration slightly overstated precision);",
            "  the prolactin update is the tighter one. The posterior reweights",
            "  a population prior built from postpartum cis literature — it does",
            "  not launder that prior's weaknesses. Use under clinical care.",
        ]
        return "\n".join(lines)


def posterior_predictive(
    result: PersonalizedResult,
    candidate_doses: List[float],
    prep_weeks: float,
    base_protocol: LactationProtocol,
    dt_weeks: float = 0.25,
) -> str:
    """Weighted supply prediction at candidate next doses — the titration
    decision table. Reuses each posterior body's own drug + gland parameters."""

    g = get_galactagogue(base_protocol.galactagogue)
    names = result.param_names
    lines = [
        "Posterior predictive — expected steady supply at candidate doses:",
        f"{'dose mg/d':>10} | {'p10':>6} | {'median':>7} | {'p90':>6}   "
        f"(mL/day, weighted over YOUR posterior)",
        "-" * 56,
    ]
    for d in candidate_doses:
        emax = result.params[:, names.index("emax")]
        ed50 = result.params[:, names.index("ed50")]
        base = result.params[:, names.index("baseline")]
        vmax = result.params[:, names.index("milk_vmax")]
        km = result.params[:, names.index("milk_km")]
        # Full re-simulation per posterior body at this dose (exact; the
        # posterior weights carry over unchanged).
        sims = np.empty(result.weights.size)
        for i in range(result.weights.size):
            m = LactationModel(
                milk_vmax_ml_day=vmax[i],
                milk_km_prolactin_ng_ml=km[i],
                maturation_tau_weeks=result.params[i, names.index("mat_tau")],
                maturation_prolactin_threshold_ng_ml=result.params[i, names.index("mat_thr")],
                ductal_tau_weeks=result.params[i, names.index("ductal_tau")],
                la_tau_weeks=result.params[i, names.index("la_tau")],
                involution_rate_per_week=result.params[i, names.index("invol")],
            )
            prl = base[i] + emax[i] * d / (ed50[i] + d)
            res = simulate_lactation(
                base_protocol, total_weeks=prep_weeks + 60, dt_weeks=dt_weeks,
                model=m, individual_prolactin_ng_ml=prl,
            )
            n_tail = max(1, int(round(8.0 / dt_weeks)))
            sims[i] = res.milk_ml_day[-n_tail:].mean()
        p10 = result.weighted_percentile(sims, 10)
        med = result.weighted_percentile(sims, 50)
        p90 = result.weighted_percentile(sims, 90)
        lines.append(f"{d:>10.1f} | {p10:>6.0f} | {med:>7.0f} | {p90:>6.0f}")
    lines.append("")
    lines.append("  Titration reading: choose the smallest dose whose posterior")
    lines.append("  median clears your supply goal — extra dose beyond the knee")
    lines.append("  buys little milk and real exposure.")
    return "\n".join(lines)


def personalize(
    agent: str,
    observations: List[Observation],
    prep_weeks: float = 30.0,
    progesterone_start_weeks: Optional[float] = None,
    pumping_sessions_per_day: float = 8.0,
    n: int = 400,
    seed: int = 7,
    spec: Optional[SensitivitySpec] = None,
    dt_weeks: float = 0.25,
) -> PersonalizedResult:
    """Reweight the population ensemble by a participant's measurements."""
    if not observations:
        raise ValueError("personalize() needs at least one Observation")
    g = get_galactagogue(agent)
    spec = spec or SensitivitySpec()
    rng = np.random.default_rng(seed)

    max_week = max(o.weeks_after_drop for o in observations)
    base_protocol = LactationProtocol(
        galactagogue=agent,
        galactagogue_dose_mg_day=float(np.mean([o.dose_mg_day for o in observations])),
        prep_weeks=prep_weeks,
        progesterone_start_weeks=progesterone_start_weeks,
        pumping_sessions_per_day=pumping_sessions_per_day,
    )

    param_names = ["emax", "ed50", "baseline", "milk_vmax", "milk_km",
                   "mat_tau", "mat_thr", "ductal_tau", "la_tau", "invol"]
    P = np.empty((n, len(param_names)))
    logw = np.empty(n)
    milk_curves = None
    weeks = None

    total_weeks = prep_weeks + max(max_week + 8.0, 52.0)

    for i in range(n):
        m = spec.sample(rng)
        emax = g.prolactin_emax_ng_ml * float(
            rng.lognormal(0.0, spec.galactagogue_cvs["prolactin_emax_ng_ml"]))
        ed50 = max(1e-3, g.prolactin_ed50_mg_day * float(
            rng.lognormal(0.0, spec.galactagogue_cvs["prolactin_ed50_mg_day"])))
        base = PROLACTIN_BASELINE_NG_ML * float(
            rng.lognormal(0.0, spec.galactagogue_cvs["baseline"]))

        P[i] = [emax, ed50, base,
                m.milk_vmax_ml_day, m.milk_km_prolactin_ng_ml,
                m.maturation_tau_weeks, m.maturation_prolactin_threshold_ng_ml,
                m.ductal_tau_weeks, m.la_tau_weeks, m.involution_rate_per_week]

        # Predict + accumulate log-likelihood per observation.
        lw = 0.0
        for o in observations:
            if o.prolactin_ng_ml is not None:
                prl = base + emax * o.dose_mg_day / (ed50 + o.dose_mg_day)
                sigma = max(0.15 * o.prolactin_ng_ml, 3.0)
                lw += -0.5 * ((prl - o.prolactin_ng_ml) / sigma) ** 2
            if o.milk_ml_day is not None:
                # Approximation: treat the observation's dose as constant from
                # the drop (documented; titration handled by re-calling with
                # per-observation doses in a future revision).
                proto = replace(base_protocol,
                                galactagogue_dose_mg_day=float(o.dose_mg_day))
                res = simulate_lactation(proto, total_weeks=total_weeks,
                                         dt_weeks=dt_weeks, model=m)
                wk_idx = int(round((prep_weeks + o.weeks_after_drop) / dt_weeks))
                pred = float(res.milk_ml_day[min(wk_idx, res.milk_ml_day.size - 1)])
                sigma = max(0.30 * o.milk_ml_day, 10.0)
                lw += -0.5 * ((pred - o.milk_ml_day) / sigma) ** 2
                if milk_curves is None:
                    milk_curves = np.empty((n, res.milk_ml_day.size))
                    weeks = res.weeks
                milk_curves[i] = res.milk_ml_day
        logw[i] = lw

    if milk_curves is None:
        # Prolactin-only observations: run the curves once at mean dose for
        # the predictive bands.
        res = simulate_lactation(base_protocol, total_weeks=total_weeks,
                                 dt_weeks=dt_weeks)
        milk_curves = np.tile(res.milk_ml_day, (n, 1))
        weeks = res.weeks

    w = np.exp(logw - logw.max())
    w /= w.sum()
    ess = 1.0 / np.sum(w ** 2) if w.sum() > 0 else 0.0

    return PersonalizedResult(
        observations=list(observations),
        weights=w,
        ess=float(ess),
        param_names=param_names,
        params=P,
        weeks=weeks,
        milk=milk_curves,
        prior_spec=spec,
    )
