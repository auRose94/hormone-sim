"""Participant-data ingest: personal inference from measured points.

This is the module that turns the model from a *population prior* into a
*personal inference tool* — the piece a titration study actually needs.

Mechanism (importance sampling, no MCMC) — the domain-generic engine lives
in :mod:`hormone_sim.measurement`; this module provides the lactation
domain's measurement model:

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
from .measurement import (
    BodySample,
    DomainPosterior,
    MeasurementDomain,
    personalize_domain,
    posterior_predictive_table,
)


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


# params row layout shared by the domain's likelihood and predictive table.
_LACTATION_PARAM_NAMES = ["emax", "ed50", "baseline", "milk_vmax", "milk_km",
                          "mat_tau", "mat_thr", "ductal_tau", "la_tau", "invol"]


def _model_from_params(params: np.ndarray) -> LactationModel:
    """Rebuild the sampled gland model from a params row (order above)."""
    return LactationModel(
        milk_vmax_ml_day=params[3],
        milk_km_prolactin_ng_ml=params[4],
        maturation_tau_weeks=params[5],
        maturation_prolactin_threshold_ng_ml=params[6],
        ductal_tau_weeks=params[7],
        la_tau_weeks=params[8],
        involution_rate_per_week=params[9],
    )


@dataclass
class LactationDomain(MeasurementDomain):
    """The lactation measurement model (see :mod:`hormone_sim.measurement`).

    Drug side: sampled prolactin Emax/ED50/baseline. Gland side: the
    ``SensitivitySpec`` sample. Measurements: prolactin (analytic Emax
    prediction) and milk volume (protocol simulation at the observation's
    dose). Also carries the population descriptors the study layer's
    shared machinery uses (predictive table, population bands/sweep).
    """

    agent: str
    base_protocol: LactationProtocol
    spec: SensitivitySpec = field(default_factory=SensitivitySpec)
    dt_weeks: float = 0.25
    total_weeks: Optional[float] = None          # personalization horizon
    include_drug_uncertainty: bool = True        # population sampling only
    threshold_ml_day: float = 100.0              # sweep criterion

    name: str = "lactation"
    param_names: List[str] = field(
        default_factory=lambda: list(_LACTATION_PARAM_NAMES))

    # predictive-table presentation (shared renderer in measurement.py)
    predictive_title: str = ("Posterior predictive — expected steady supply "
                             "at candidate doses:")
    predictive_candidate_header: str = f"{'dose mg/d':>10}"
    predictive_metrics: List[Tuple[str, int]] = field(default_factory=lambda: [
        ("p10", 6), ("median", 7), ("p90", 6)])
    predictive_suffix: str = "(mL/day, weighted over YOUR posterior)"
    predictive_rule: str = "-" * 56
    predictive_notes: Tuple[str, ...] = (
        "  Titration reading: choose the smallest dose whose posterior",
        "  median clears your supply goal — extra dose beyond the knee",
        "  buys little milk and real exposure.",
    )

    def __post_init__(self) -> None:
        self.g = get_galactagogue(self.agent)

    # --- personalization -----------------------------------------------------
    def sample_body(self, rng: np.random.Generator,
                    observations) -> BodySample:
        m = self.spec.sample(rng)
        emax = self.g.prolactin_emax_ng_ml * float(
            rng.lognormal(0.0, self.spec.galactagogue_cvs["prolactin_emax_ng_ml"]))
        ed50 = max(1e-3, self.g.prolactin_ed50_mg_day * float(
            rng.lognormal(0.0, self.spec.galactagogue_cvs["prolactin_ed50_mg_day"])))
        base = PROLACTIN_BASELINE_NG_ML * float(
            rng.lognormal(0.0, self.spec.galactagogue_cvs["baseline"]))
        params = np.array([emax, ed50, base,
                           m.milk_vmax_ml_day, m.milk_km_prolactin_ng_ml,
                           m.maturation_tau_weeks,
                           m.maturation_prolactin_threshold_ng_ml,
                           m.ductal_tau_weeks, m.la_tau_weeks,
                           m.involution_rate_per_week])
        return BodySample(params=params)

    def log_likelihood(self, body: BodySample, obs) -> float:
        emax, ed50, base = body.params[0], body.params[1], body.params[2]
        lw = 0.0
        if obs.prolactin_ng_ml is not None:
            prl = base + emax * obs.dose_mg_day / (ed50 + obs.dose_mg_day)
            sigma = max(0.15 * obs.prolactin_ng_ml, 3.0)
            lw += -0.5 * ((prl - obs.prolactin_ng_ml) / sigma) ** 2
        if obs.milk_ml_day is not None:
            # Approximation: treat the observation's dose as constant from
            # the drop (documented; titration handled by re-calling with
            # per-observation doses in a future revision).
            proto = replace(self.base_protocol,
                            galactagogue_dose_mg_day=float(obs.dose_mg_day))
            res = simulate_lactation(proto, total_weeks=self.total_weeks,
                                     dt_weeks=self.dt_weeks,
                                     model=_model_from_params(body.params))
            wk_idx = int(round((self.base_protocol.prep_weeks
                                + obs.weeks_after_drop) / self.dt_weeks))
            pred = float(res.milk_ml_day[min(wk_idx, res.milk_ml_day.size - 1)])
            sigma = max(0.30 * obs.milk_ml_day, 10.0)
            lw += -0.5 * ((pred - obs.milk_ml_day) / sigma) ** 2
            body.curves["milk"] = res.milk_ml_day
            body.weeks = res.weeks
        return lw

    def finalize(self, bodies: List[BodySample], observations) -> None:
        """Prolactin-only logs: fill every body's milk curve with the
        population-default simulation at the mean dose (the predictive
        bands then carry the prior, weightable by the prolactin update)."""
        if all("milk" in b.curves for b in bodies):
            return
        res = simulate_lactation(self.base_protocol,
                                 total_weeks=self.total_weeks,
                                 dt_weeks=self.dt_weeks)
        for b in bodies:
            if "milk" not in b.curves:
                b.curves["milk"] = res.milk_ml_day
                b.weeks = res.weeks

    def curve_names(self) -> Tuple[str, ...]:
        return ("milk",)

    # --- posterior-predictive table ------------------------------------------
    def candidate_label(self, candidate: float) -> str:
        return f"{candidate:>10.1f}"

    def predict_candidate(self, params: np.ndarray,
                          candidate: float) -> List[float]:
        m = _model_from_params(params)
        emax, ed50, base = params[0], params[1], params[2]
        prl = base + emax * candidate / (ed50 + candidate)
        res = simulate_lactation(
            self.base_protocol, total_weeks=self.base_protocol.prep_weeks + 60,
            dt_weeks=self.dt_weeks, model=m, individual_prolactin_ng_ml=prl)
        n_tail = max(1, int(round(8.0 / self.dt_weeks)))
        return [float(res.milk_ml_day[-n_tail:].mean())]

    # --- population bands / sweep ---------------------------------------------
    def population_body(self, rng: np.random.Generator,
                        candidate: float) -> BodySample:
        m = self.spec.sample(rng)
        p_drug = None
        if self.include_drug_uncertainty:
            emax = self.g.prolactin_emax_ng_ml * float(
                rng.lognormal(0.0, self.spec.galactagogue_cvs["prolactin_emax_ng_ml"]))
            ed50 = max(1e-3, self.g.prolactin_ed50_mg_day * float(
                rng.lognormal(0.0, self.spec.galactagogue_cvs["prolactin_ed50_mg_day"])))
            base = PROLACTIN_BASELINE_NG_ML * float(
                rng.lognormal(0.0, self.spec.galactagogue_cvs["baseline"]))
            d = max(0.0, float(candidate))
            p_drug = base + emax * d / (ed50 + d)
        proto = replace(self.base_protocol,
                        galactagogue_dose_mg_day=float(candidate))
        res = simulate_lactation(proto, total_weeks=self.total_weeks,
                                 dt_weeks=self.dt_weeks, model=m,
                                 individual_prolactin_ng_ml=p_drug)
        params = np.array([0.0, 0.0, 0.0, m.milk_vmax_ml_day,
                           m.milk_km_prolactin_ng_ml, m.maturation_tau_weeks,
                           m.maturation_prolactin_threshold_ng_ml,
                           m.ductal_tau_weeks, m.la_tau_weeks,
                           m.involution_rate_per_week])
        body = BodySample(params=params, weeks=res.weeks)
        body.curves["milk"] = res.milk_ml_day
        return body

    def sweep_summary(self, curves: dict) -> dict:
        milk = curves["milk"]
        n_tail = max(1, int(round(8.0 / self.dt_weeks)))
        return {"steady_milk": float(milk[-n_tail:].mean())}

    def sweep_hit(self, summary: dict) -> bool:
        return summary["steady_milk"] >= self.threshold_ml_day


@dataclass
class PersonalizedResult(DomainPosterior):
    """Posterior over modeled bodies given a participant's observations.

    The lactation view of :class:`~hormone_sim.measurement.DomainPosterior`
    (same weighting machinery, milk-shaped accessors).
    """

    prior_spec: SensitivitySpec

    @property
    def milk(self) -> np.ndarray:
        return self.curves["milk"]

    def band(self, q: float = 0.8) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Weighted (lower, median, upper) milk curves for central mass q."""
        return super().band(q, which="milk")

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
    decision table. Reuses each posterior body's own drug + gland parameters;
    rendered by the domain-generic table machinery
    (:func:`hormone_sim.measurement.posterior_predictive_table`)."""
    domain = LactationDomain(
        agent=base_protocol.galactagogue, base_protocol=base_protocol,
        dt_weeks=dt_weeks,
    )
    return posterior_predictive_table(domain, result, candidate_doses)


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
    spec = spec or SensitivitySpec()

    max_week = max(o.weeks_after_drop for o in observations)
    base_protocol = LactationProtocol(
        galactagogue=agent,
        galactagogue_dose_mg_day=float(np.mean([o.dose_mg_day for o in observations])),
        prep_weeks=prep_weeks,
        progesterone_start_weeks=progesterone_start_weeks,
        pumping_sessions_per_day=pumping_sessions_per_day,
    )
    total_weeks = prep_weeks + max(max_week + 8.0, 52.0)

    domain = LactationDomain(
        agent=agent, base_protocol=base_protocol, spec=spec,
        dt_weeks=dt_weeks, total_weeks=total_weeks,
    )
    post = personalize_domain(domain, observations, n=n, seed=seed)
    return PersonalizedResult(
        domain=domain, observations=post.observations, weights=post.weights,
        ess=post.ess, param_names=post.param_names, params=post.params,
        curves=post.curves, weeks=post.weeks, prior_spec=spec,
    )
