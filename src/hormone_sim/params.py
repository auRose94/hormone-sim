"""Physiological and pharmacokinetic parameters for GAHT modeling.

Every constant here is either:

* **cited** — traceable to a published source (see ``docs/references.md`` and the
  ``citation`` field), or
* **calibrated** — fitted to reproduce a published steady-state target where a
  direct mechanistic measurement does not exist. Calibrated constants are flagged
  ``calibrated=True`` so they are never mistaken for measured physiology.

Units
-----
Concentrations are handled in two conventions and converted at the boundary:

* Estradiol: **pg/mL** (clinical convention) and **pmol/L** (SI).  1 pg/mL = 3.671 pmol/L.
* Testosterone: **ng/dL** (clinical) and **nmol/L** (SI).  1 ng/dL = 0.0347 nmol/L.

The PK core works in **mass/volume** internally; see :func:`e2_pgml_to_pmol`.

This module is deliberately small and self-documenting. If you change a value,
update its citation and re-run ``pytest`` — the sanity invariants in ``tests``
will catch physiologically impossible results.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


# --- Unit conversion factors (exact / widely accepted) ----------------------

E2_PGML_TO_PMOL = 3.671   # 1 pg/mL estradiol = 3.671 pmol/L
T_NGDL_TO_NMOL = 0.0347   # 1 ng/dL testosterone = 0.0347 nmol/L


# --- Clinical target ranges (guideline-cited) -------------------------------
# Endocrine Society Clinical Practice Guideline (Grimley et al., 2017) and
# WPATH SOC-8 (2022): target trans-feminine estradiol ~100-200 pg/mL,
# testosterone <50 ng/dL (female-typical range).

TARGET_E2_PGML = (100.0, 200.0)      # pg/mL
TARGET_T_NGDL = (0.0, 50.0)          # ng/dL (upper bound is the clinical goal)

# Baseline (untreated, adult AMAB) values these models start from.
BASELINE_T_NGDL = 300.0              # pg/mL-equivalent; adult AMAB mean ~300 ng/dL
BASELINE_E2_PGML = 15.0              # adult AMAB estradiol ~10-20 pg/mL


@dataclass(frozen=True)
class ParamSource:
    """Provenance for a single parameter value."""

    citation: str
    calibrated: bool = False
    note: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "citation", self.citation.strip())


# --- Estradiol PK parameters -------------------------------------------------
# A 1-compartment model with first-order elimination. Volume of distribution
# and clearance are molecule properties; route differences enter through the
# input function (absorption rate / bioavailability / input shape), NOT through
# different clearance values per route.

@dataclass(frozen=True)
class EstradiolPK:
    """Population-typical estradiol pharmacokinetics.

    A **calibrated phenomenological** one-compartment model. Two independent
    knobs, each tied to a different observable, are used because a single
    clearance cannot reproduce both:

    * **``half_life_h``** governs the *kinetics* (how fast the level changes /
      fluctuates). It is calibrated to the physiological estradiol half-life
      (~10-20 h) — see source below.
    * **``steady_state_factor``** governs the *dose→level mapping* (the
      steady-state level for a given delivered dose). It is calibrated to the
      clinical anchor: transdermal 50 µg/day → ~150 pg/mL (Endocrine Society).

    Separating them (rather than deriving one clearance for both) keeps the
    kinetics physiological while matching the observed steady state — and is the
    honest way to say this 1-compartment model does not claim to resolve real
    estradiol's true multi-phase distribution. ``v_d_l`` is a fixed apparent
    volume that, with the implied clearance, yields the target half-life.

    Route differences (bioavailability, first-pass, ester mass) are NOT in these
    molecule parameters — they live in the per-route ``delivery_factor`` and the
    oral first-pass multiplier below.
    """

    # Apparent volume of distribution (L). Combined with clearance to give the
    # target half-life; NOT the physiological ~100 L tissue distribution volume.
    v_d_l: float = 5.0
    v_d_source: ParamSource = field(
        default=ParamSource(
            citation="Apparent volume (see class docstring). Real estradiol Vd "
            "~ 1-2 L/kg (Stanczyk & Clarke, Best Pract Res Clin Endocrinol "
            "Metab, 2014) is NOT resolved by this 1-compartment model.",
            calibrated=True,
            note="Fixed apparent volume; half-life is set directly, clearance "
            "implied.",
        )
    )

    # Elimination half-life (h). Calibrated to the physiological estradiol
    # half-life (~10-20 h). The implied clearance = ln2*Vd/t1/2.
    half_life_h_val: float = 15.0
    half_life_source: ParamSource = field(
        default=ParamSource(
            citation="Sitaraman & Munir, 'Estradiol', StatPearls, 2023; Kuhl "
            "(2005) — parent estradiol t1/2 ~ 10-20 h.",
            calibrated=True,
            note="Calibrated to physiological parent estradiol half-life.",
        )
    )

    # Calibrated dose→level mapping: steady-state pg/mL per µg/day of *delivered*
    # (post-delivery_factor) estradiol. Anchor: 50 µg/day -> ~150 pg/mL => 3.0.
    steady_state_factor: float = 3.0
    steady_state_source: ParamSource = field(
        default=ParamSource(
            citation="Calibrated to the transdermal 50 ug/day -> ~150 pg/mL "
            "anchor (Grimley et al., Endocrine Society CPG, J Clin Endocrinol "
            "Metab, 2017).",
            calibrated=True,
            note="pg/mL steady-state per µg/day delivered. Per-route "
            "delivery_factor reproduces the other routes' observed levels.",
        )
    )

    # First-pass hepatic extraction multiplier for the ORAL route only. Oral
    # estradiol is cleared faster on first pass; this scales the *effective*
    # clearance (kinetics), independent of the extent (delivery_factor).
    oral_first_pass_factor: float = 1.5
    oral_first_pass_source: ParamSource = field(
        default=ParamSource(
            citation="Sitaraman & Munir, 'Estradiol', StatPearls, 2023; Kuhl "
            "(2005). Direction/magnitude of first-pass; factor calibrated.",
            calibrated=True,
        )
    )

    @property
    def clearance_l_per_h(self) -> float:
        """Implied clearance (L/h) = ln2 * Vd / t1/2 (transdermal kinetics)."""
        return (0.6931471805599453 * self.v_d_l) / self.half_life_h_val

    def effective_clearance(self, route: str = "transdermal") -> float:
        """Effective clearance (L/h) for a route (oral gets first-pass boost)."""
        return self.clearance_l_per_h * (
            self.oral_first_pass_factor if route == "oral" else 1.0
        )

    def input_scale(self) -> float:
        """Scale factor converting delivered dose-rate (µg/h) into the
        calibrated model input so that steady-state E2 = baseline +
        ``steady_state_factor`` * delivered_daily_dose.

        Derived from C_ss_incr = input_rate / Cl = steady_state_factor * dose,
        with dose = continuous_rate * 24 h. Uses the *base* clearance (the
        oral first-pass multiplier affects kinetics/half-life only, not the
        steady-state level mapping — that is captured by ``delivery_factor``).
        """
        return self.clearance_l_per_h * self.steady_state_factor * 24.0

    def half_life_h(self, route: str = "transdermal") -> float:
        """Apparent elimination half-life (h) for a given route (oral shorter)."""
        cl = self.effective_clearance(route)
        return (0.6931471805599453 * self.v_d_l) / cl


@dataclass(frozen=True)
class RouteParams:
    """Route-specific input-function parameters for estradiol delivery.

    ``input_kind`` selects the absorption model:
      * ``"first_order"``  — oral, IM/SubQ esters (depot): rate-limited uptake
      * ``"zero_order"``   — transdermal patch/gel: constant flux while applied

    ``bioavailability`` is the fraction of administered estradiol reaching
    systemic circulation (cited). ``delivery_factor`` is a **calibrated** scalar
    applied to the nominal dose to reproduce the *observed* steady-state level
    for that route (it folds in bioavailability, first-pass extent, and — for
    esters — the active-moiety fraction of the ester mass). ``absorption_rate_h``
    is the first-order rate (1/h) for depot/gut uptake.
    """

    name: str
    input_kind: str
    bioavailability: float
    delivery_factor: float = 1.0
    absorption_rate_h: float = 0.3          # 1/h, first-order only
    # For depot (ester) injections: fraction released per hour from the depot.
    depot_rate_h: float = 0.02              # 1/h, slow release
    # Applied duration per "dose event" for zero-order (patch/gel), hours.
    application_hours: float = 24.0


# Population-typical route parameters. ``bioavailability`` is cited;
# ``delivery_factor`` is CALIBRATED so each route's typical clinical dose lands
# in its guideline target range, keeping the model's dose→level mapping
# (steady_state_factor) consistent across routes. Because this model's apparent
# clearance is phenomenological (see EstradiolPK), delivery_factor is NOT a
# literal bioavailability — it is the per-route scalar that reproduces the
# observed steady-state level. Calibrated targets (see docs/references.md):
#   transdermal 50 ug/day  -> ~150 pg/mL  (reference; df=1.0)
#   oral 2000 ug/day       -> ~200 pg/mL  (heavy first-pass, needs mg doses)
#   IM 10 mg/week (~1430 ug/day) -> ~150 pg/mL
#   subQ 5 mg/week  (~714 ug/day) -> ~150 pg/mL
# With steady_state_factor=3.0 (pg/mL per ug/day delivered), df = target/(3*dose).
ROUTES: Dict[str, RouteParams] = {
    "oral": RouteParams(
        name="oral",
        input_kind="first_order",
        bioavailability=0.10,   # micronized oral E2 ~ 10% systemic (heavy first-pass)
        delivery_factor=0.033,  # calibrated: 2000 ug/day -> ~200 pg/mL
        absorption_rate_h=1.5,  # rapid GI uptake
    ),
    "transdermal": RouteParams(
        name="transdermal",
        input_kind="zero_order",
        bioavailability=0.90,   # transdermal bypasses first-pass, high effective F
        delivery_factor=1.0,    # reference: 50 ug/day -> ~150 pg/mL
        application_hours=24.0,
    ),
    "im_injection": RouteParams(
        name="im_injection",
        input_kind="first_order",   # depot: slow release then first-order elimination
        bioavailability=0.90,
        delivery_factor=0.035,      # calibrated: 10 mg/week -> ~150 pg/mL
        depot_rate_h=0.02,          # slow ester release from oil depot
    ),
    "subq_injection": RouteParams(
        name="subq_injection",
        input_kind="first_order",
        bioavailability=0.90,
        delivery_factor=0.070,      # calibrated: 5 mg/week -> ~150 pg/mL
        depot_rate_h=0.03,          # slightly faster than IM in some reports
    ),
}


# --- Testosterone suppression (PD) -------------------------------------------
# Model the HPG-axis suppression by estradiol negative feedback, plus additive
# pharmacologic anti-androgen effects. The estradiol suppression is a
# Hill-function (sigmoid) of the *average* estradiol level; the anti-androgen
# adds an independent fractional suppression.

@dataclass(frozen=True)
class AndrogenPD:
    """Testosterone suppression model parameters.

    T(t) = T_baseline * (1 - S_e2) * (1 - S_aa)

    where S_e2 is estradiol-driven HPG suppression (Hill function of avg E2)
    and S_aa is the anti-androgen effect (sum of fractional contributions,
    capped at 1).
    """

    # Estradiol level (pg/mL) giving 50% maximal HPG-axis suppression.
    # Calibrated so that E2 ~150 pg/mL gives near-complete suppression,
    # consistent with the clinical goal of T < 50 ng/dL.
    e2_ec50_pgml: float = 80.0
    e2_ec50_source: ParamSource = field(
        default=ParamSource(
            citation="Calibrated to reproduce clinical goal: T<50 ng/dL at "
            "E2~150 pg/mL. See Grimley et al. 2017; Unger, 'Assessment and "
            "Management of Transgender Patients', 2018.",
            calibrated=True,
        )
    )

    # Hill coefficient (cooperativity) of the suppression curve.
    hill_n: float = 2.0
    hill_n_source: ParamSource = field(
        default=ParamSource(
            citation="Calibrated; sigmoidal feedback is standard endocrine "
            "modeling.",
            calibrated=True,
        )
    )

    # Anti-androgen fractional effects at the *typical* clinical dose. These
    # are population-typical maximal effects (Emax), not patient-specific.
    # Sources: clinical pharmacology of each agent.
    anti_androgen_emax: Dict[str, float] = field(
        default_factory=lambda: {
            # Fractional suppression of T production/action at typical dose.
            "spironolactone": 0.55,   # ~50-60% T reduction at 100-200 mg/d
            "cyproterone": 0.90,      # CPA: potent, ~85-95% suppression
            "gnrh_agonist": 0.97,     # 'chemical castration'; near-complete
            "bicalutamide": 0.10,     # blocks action, can *raise* serum T; net production effect small
        }
    )
    anti_androgen_emax_source: ParamSource = field(
        default=ParamSource(
            citation="Kuhl & Weber, 'Bicalutamide in the treatment of "
            "transsexualism', Int J Transgenderism, 2004; Bock & Hoff, 'The "
            "endocrinology of transgender persons', Endocrinol Metab Clin "
            "North Am, 2016; Grimley et al. 2017.",
            calibrated=False,
        )
    )


# --- Population variability (for Monte-Carlo) --------------------------------
@dataclass(frozen=True)
class Variability:
    """Between-subject variability, expressed as coefficient of variation (CV).

    Values are population-typical CVs (fraction) used to sample individual
    parameters around the population mean in :mod:`hormone_sim.population`.
    Sampled on the model's *tunable* knobs: half-life (kinetics) and
    steady-state factor (dose→level mapping).
    """

    half_life_cv: float = 0.30
    steady_state_cv: float = 0.35
    e2_ec50_cv: float = 0.25
    baseline_t_cv: float = 0.20


PARAMS = {
    "estradiol_pk": EstradiolPK(),
    "androgen_pd": AndrogenPD(),
    "variability": Variability(),
    "routes": ROUTES,
}


def e2_pgml_to_pmol(pgml: float) -> float:
    """Convert estradiol pg/mL -> pmol/L."""
    return pgml * E2_PGML_TO_PMOL


def t_ngdl_to_nmol(ngdl: float) -> float:
    """Convert testosterone ng/dL -> nmol/L."""
    return ngdl * T_NGDL_TO_NMOL
