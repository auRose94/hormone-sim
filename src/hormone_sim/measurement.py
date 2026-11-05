"""Measurement-model abstraction: the study layer's domain-generic engine.

The study layer used to be lactation-shaped all the way down: the
importance-sampling personalization loop, the posterior-predictive table,
and the population sweep each carried milk-and-prolactin assumptions in
their code bodies. Adding the trans-masculine domain meant copying the
machinery and renaming the outcomes — which is how silent divergence
happens. This module is the fix: the *machinery* lives here once, and each
domain (lactation, T logs, ...) describes itself through a small adapter —
the **measurement model**:

* :meth:`MeasurementDomain.sample_body` — draw one modeled body's parameters
  (and any curves it can pre-compute) from the population prior;
* :meth:`MeasurementDomain.log_likelihood` — the measurement model proper:
  how likely one observed datapoint is under that body (each measurement
  kind — prolactin lab, milk diary, T draw, CBC — gets a pragmatic sigma);
* :meth:`MeasurementDomain.curve_names` — which per-body outcome curves the
  posterior keeps;
* optional ``finalize(bodies, observations)`` — a post-pass (e.g. the T
  domain's dominant-regimen Hct curves);
* predictive-table descriptors (title/columns/notes + ``predict_candidate``)
  and population-sweep descriptors (``population_body``/``sweep_summary``/
  ``sweep_hit``) — see the respective functions below.

What the engine guarantees, for every domain equally:

* importance sampling with the exact weighting math (softmax-stabilized
  log-weights, normalized; ESS = 1/sum(w^2));
* bands over points everywhere (weighted percentiles of the ensemble ARE
  the posterior; unweighted percentiles are the population prior);
* the ESS diagnostic honestly stating when posterior ≈ prior.

Determinism contract: ``personalize_domain`` draws RNG numbers ONLY through
``sample_body`` (in body order), then evaluates likelihoods. Domains must
never draw inside ``log_likelihood`` — the simulators are deterministic, so
this keeps any domain's results reproducible and lets the per-public-API
wrappers reproduce their pre-abstraction behavior bit-for-bit.

Not medical advice: the engine is analysis machinery for clinician-
supervised titration studies; the per-domain caveats live with the domains
and must be carried into every derivative.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


# --- Domain interface ------------------------------------------------------------

class MeasurementDomain:
    """Structural interface for a study-layer domain (see module docstring).

    Concrete domains (e.g. :class:`~hormone_sim.personalization.LactationDomain`,
    :class:`~hormone_sim.personalization_androgen.AndrogenDomain`) implement
    the required members; they do not need to inherit this class. Optional
    hooks (``finalize``, the predictive-table and population-sweep
    descriptors) are discovered with ``getattr`` so a domain implements only
    what it supports.
    """

    name: str
    param_names: List[str]

    def sample_body(self, rng: np.random.Generator,
                    observations: Sequence[Any]) -> "BodySample":
        raise NotImplementedError

    def log_likelihood(self, body: "BodySample", obs: Any) -> float:
        raise NotImplementedError

    def curve_names(self) -> Tuple[str, ...]:
        raise NotImplementedError


@dataclass
class BodySample:
    """One modeled body: sampled parameters + cached predictions/curves.

    ``params`` is a row aligned with the domain's ``param_names``. ``curves``
    holds per-body outcome curves over ``weeks`` (filled either at sampling
    time, inside likelihood evaluation, or by ``finalize`` — the engine only
    reads them after all three have run).
    """

    params: np.ndarray
    curves: Dict[str, np.ndarray] = field(default_factory=dict)
    weeks: Optional[np.ndarray] = None


# --- Posterior -------------------------------------------------------------------

@dataclass
class DomainPosterior:
    """Importance-sampling posterior over modeled bodies, domain-generic.

    The weighted percentiles of the ensemble ARE the posterior; the
    unweighted percentiles are the population prior (a like-for-like
    comparison — same bodies, weights removed). ESS diagnostics say whether
    the observations actually moved the ensemble.
    """

    domain: Any
    observations: List[Any]
    weights: np.ndarray            # normalized importance weights, shape (n,)
    ess: float
    param_names: List[str]
    params: np.ndarray             # (n, k) sampled parameter values
    curves: Dict[str, np.ndarray]  # name -> (n, t) unweighted outcome curves
    weeks: Optional[np.ndarray]

    def weighted_percentile(self, values: np.ndarray, p: float) -> float:
        idx = np.argsort(values)
        v = values[idx]
        w = self.weights[idx]
        cdf = np.cumsum(w)
        cdf /= cdf[-1]
        return float(np.interp(p / 100.0, cdf, v))

    def band(self, q: float = 0.8, which: str = None
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Weighted (lower, median, upper) curves for central mass ``q``."""
        if which is None:
            which = next(iter(self.curves))
        matrix = self.curves[which]
        lo_pct = (1.0 - q) / 2.0 * 100.0
        hi_pct = (1.0 + q) / 2.0 * 100.0
        n_t = matrix.shape[1]
        lo = np.array([self.weighted_percentile(matrix[:, j], lo_pct)
                       for j in range(n_t)])
        med = np.array([self.weighted_percentile(matrix[:, j], 50.0)
                        for j in range(n_t)])
        hi = np.array([self.weighted_percentile(matrix[:, j], hi_pct)
                       for j in range(n_t)])
        return lo, med, hi

    def posterior_param(self, name: str, p: float = 50.0) -> float:
        return self.weighted_percentile(
            self.params[:, self.param_names.index(name)], p)


def personalize_domain(domain: MeasurementDomain,
                       observations: Sequence[Any],
                       n: int, seed: int) -> DomainPosterior:
    """Reweight a domain's population ensemble by a participant's observations.

    The domain-generic importance-sampling engine. RNG draws happen only in
    ``sample_body`` (determinism contract, module docstring); likelihoods
    and any ``finalize`` post-pass are pure computation.
    """
    if not observations:
        raise ValueError(f"{getattr(domain, 'name', type(domain).__name__)} "
                         "personalization needs at least one observation")
    rng = np.random.default_rng(seed)

    bodies = [domain.sample_body(rng, observations) for _ in range(n)]

    logw = np.empty(n)
    for i, body in enumerate(bodies):
        lw = 0.0
        for obs in observations:
            lw += domain.log_likelihood(body, obs)
        logw[i] = lw

    finalize = getattr(domain, "finalize", None)
    if finalize is not None:
        finalize(bodies, observations)

    curves: Dict[str, np.ndarray] = {}
    for name in domain.curve_names():
        filled = [b.curves.get(name) for b in bodies]
        if any(c is None for c in filled):
            missing = sum(1 for c in filled if c is None)
            raise ValueError(
                f"domain {getattr(domain, 'name', type(domain).__name__)!r}: "
                f"{missing}/{n} bodies are missing curve {name!r} after "
                "sampling/likelihood/finalize — the domain must fill every "
                "curve for every body")
        curves[name] = np.stack(filled)

    weeks = next((b.weeks for b in bodies if b.weeks is not None), None)

    w = np.exp(logw - logw.max())
    w /= w.sum()
    ess = 1.0 / np.sum(w ** 2) if w.sum() > 0 else 0.0

    return DomainPosterior(
        domain=domain,
        observations=list(observations),
        weights=w,
        ess=float(ess),
        param_names=list(domain.param_names),
        params=np.vstack([b.params for b in bodies]),
        curves=curves,
        weeks=weeks,
    )


# --- Posterior-predictive table ---------------------------------------------------
# Domain descriptors used here:
#   predictive_title: str                       first line
#   predictive_candidate_header: str            left-hand column(s) of the header
#   predictive_metrics: List[Tuple[str, int]]   (column name, width) per metric
#   predictive_suffix: str                      units note after the header
#   predictive_rule: str                        the "-" divider line
#   candidate_label(candidate) -> str           left-hand cell(s) of a row
#   predict_candidate(params_row, candidate)    per-body metric values
#   predictive_notes: Sequence[str]             the "titration reading" lines

def posterior_predictive_table(domain: MeasurementDomain,
                               posterior: DomainPosterior,
                               candidates: Sequence[Any]) -> str:
    """Weighted posterior-median predictions at candidate next regimens —
    the titration decision table, rendered uniformly for every domain.

    Each candidate is re-simulated per posterior body (the weights carry
    over unchanged), and the table shows weighted medians. The columns are
    the domain's choice (milk percentiles vs trough/mean/peak + Hct
    plateau); the machinery — weighting, formatting, the reading guidance —
    is shared.
    """
    lines = [
        domain.predictive_title,
        " | ".join([domain.predictive_candidate_header]
                   + [f"{name:>{w}}" for name, w in domain.predictive_metrics])
        + (f"   {domain.predictive_suffix}" if domain.predictive_suffix else ""),
        domain.predictive_rule,
    ]
    for candidate in candidates:
        vals = np.empty((posterior.weights.size,
                         len(domain.predictive_metrics)))
        for i in range(posterior.weights.size):
            vals[i] = domain.predict_candidate(posterior.params[i], candidate)
        cells = [f"{posterior.weighted_percentile(vals[:, m], 50):>{w}.0f}"
                 for m, (_name, w) in enumerate(domain.predictive_metrics)]
        lines.append(f"{domain.candidate_label(candidate)} | "
                     + " | ".join(cells))
    lines.append("")
    lines.extend(domain.predictive_notes)
    return "\n".join(lines)


# --- Population bands & sweep ------------------------------------------------------

def population_bands(domain: MeasurementDomain, candidate: Any,
                     n: int, seed: int) -> Tuple[Optional[np.ndarray],
                                                 Dict[str, np.ndarray]]:
    """Monte-Carlo a domain's population at one candidate regimen.

    Requires the optional ``population_body(rng, candidate)`` hook. Returns
    (weeks, {curve_name: (n, t) array}). The same paired-seed convention as
    the sweep: callers use one seed per candidate so the same "bodies" are
    evaluated across candidates.
    """
    population_body = getattr(domain, "population_body", None)
    if population_body is None:
        raise TypeError(f"{type(domain).__name__} does not implement "
                        "population_body; it cannot produce population bands")
    rng = np.random.default_rng(seed)
    bodies = [population_body(rng, candidate) for _ in range(n)]
    weeks = next((b.weeks for b in bodies if b.weeks is not None), None)
    curves = {name: np.stack([b.curves[name] for b in bodies])
              for name in domain.curve_names()}
    return weeks, curves


@dataclass
class DomainSweepRow:
    """One candidate's population outcome, between-body heterogeneity explicit.

    ``metrics`` maps each of the domain's summary metrics to its
    (p10, p50, p90) across modeled bodies; ``frac_hit`` is the fraction of
    bodies meeting the domain's criterion. Bands over points, as everywhere
    in this package.
    """

    label: str
    metrics: Dict[str, Tuple[float, float, float]]
    frac_hit: float

    def metric(self, name: str, p: str = "p50") -> float:
        idx = {"p10": 0, "p50": 1, "p90": 2}[p]
        return self.metrics[name][idx]


def population_sweep(domain: MeasurementDomain,
                     candidates: Sequence[Any],
                     n: int, seed: int) -> List[DomainSweepRow]:
    """Sweep candidate regimens, quantifying between-body heterogeneity.

    The same seed is used for every candidate, so the same modeled bodies
    are evaluated at each — rows are paired, and population medians move
    smoothly across the sweep instead of re-shuffling. Requires the
    ``population_body`` / ``sweep_summary`` / ``sweep_hit`` hooks.
    """
    for hook in ("population_body", "sweep_summary", "sweep_hit"):
        if getattr(domain, hook, None) is None:
            raise TypeError(f"{type(domain).__name__} does not implement "
                            f"{hook}; it cannot be swept")
    rows: List[DomainSweepRow] = []
    for candidate in candidates:
        rng = np.random.default_rng(seed)   # paired samples across candidates
        summaries: List[Dict[str, float]] = []
        hits: List[bool] = []
        for _ in range(n):
            body = domain.population_body(rng, candidate)
            s = domain.sweep_summary(body.curves)
            summaries.append(s)
            hits.append(bool(domain.sweep_hit(s)))
        names = list(summaries[0].keys())
        metrics: Dict[str, Tuple[float, float, float]] = {}
        for name in names:
            vals = np.array([s[name] for s in summaries])
            metrics[name] = (float(np.percentile(vals, 10)),
                             float(np.percentile(vals, 50)),
                             float(np.percentile(vals, 90)))
        rows.append(DomainSweepRow(
            label=domain.candidate_label(candidate),
            metrics=metrics,
            frac_hit=float(np.mean(hits)),
        ))
    return rows
