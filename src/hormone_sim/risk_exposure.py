"""Risk-exposure framing for D2-antagonist galactagogue use.

This module does NOT output a per-person probability of tardive dyskinesia —
no honest model can, because per-person incidence for short-course,
low-dose, monitored use in young adults is simply not established. What it
does do is structure the exposure the way the actual risk factors are known
to work, so a clinician and patient can see what is and isn't being assumed.

**The risk structure (what is known):**

* Tardive dyskinesia (TD) risk rises with **duration** and **cumulative dose**.
  FDA labeling for metoclopramide: use >12 weeks is not recommended outside
  rare cases — duration is the dominant, modifiable factor.
* **Age is the strongest single risk factor.** Incidence estimates from
  chronic-use populations (which skew elderly) do not transfer to a young,
  short-course patient. This is genuine risk heterogeneity, not reassurance.
* **No biomarker exists.** Unlike thrombotic risk (CBC/coag monitoring), TD
  monitoring is a structured clinical exam — the AIMS (Abnormal Involuntary
  Movement Scale). TD can also emerge or worsen *after* the drug is stopped,
  so monitoring must extend past the exposure window.
* Rare short-course cases are documented; the floor is not zero.

**The asymmetry vs thrombotic risk (the honest version of "is it melodramatic"):**

Estradiol's VTE risk is acute, observable on labs, reversible on stopping, and
strongly route/dose-dependent — which is why a monitored patient on
transdermal estradiol can defensibly treat population warnings as
not-applying-to-them. TD risk is cumulative, observable only by clinical
exam, potentially irreversible, and duration/cumulative-dose-dependent. The
*mitigation logic* is the same in both cases — lowest effective exposure,
shortest effective duration, structured monitoring — but the failure mode is
less forgiving. Both facts are true at once.

Tiers here are **framing, not prophecy**: "low" means *the known modifiers
are all in the favorable direction*, not *risk is zero*.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

from .galactagogues import get_galactagogue


@dataclass(frozen=True)
class ExposureProfile:
    """Structured exposure assessment for one galactagogue course."""

    agent: str
    dose_mg_day: float
    weeks: float
    age: int
    cumulative_mg: float
    duration_tier: str
    age_tier: str
    cumulative_tier: str
    overall_tier: str
    factors: List[str] = field(default_factory=list)
    monitoring: List[str] = field(default_factory=list)
    notes: str = ""

    def report(self) -> str:
        lines = [
            f"Exposure profile — {self.agent} {self.dose_mg_day:g} mg/day "
            f"for {self.weeks:g} weeks, age {self.age}",
            f"  Cumulative exposure: {self.cumulative_mg:,.0f} mg "
            f"({self.cumulative_tier})",
            f"  Duration tier:   {self.duration_tier}",
            f"  Age tier:        {self.age_tier}",
            f"  Overall framing: {self.overall_tier.upper()}",
        ]
        if self.factors:
            lines.append("  Known modifiers present:")
            for f in self.factors:
                lines.append(f"    - {f}")
        lines.append("  Monitoring plan:")
        for m in self.monitoring:
            lines.append(f"    - {m}")
        if self.notes:
            lines.append(f"  Notes: {self.notes}")
        return "\n".join(lines)


def _duration_tier(weeks: float) -> str:
    if weeks <= 4:
        return "low (short course)"
    if weeks <= 12:
        return "moderate (approaching the FDA 12-week cap)"
    return "HIGH (exceeds the FDA 12-week recommendation)"


def _age_tier(age: int) -> str:
    if age < 40:
        return "low (strongest TD risk factor is advanced age; not present)"
    if age < 65:
        return "elevated"
    return "HIGH (age is the dominant TD risk factor)"


def _cumulative_tier(cumulative_mg: float) -> str:
    if cumulative_mg < 2000:
        return "low"
    if cumulative_mg <= 6000:
        return "moderate"
    return "elevated"


def exposure_report(
    agent: str = "metoclopramide",
    dose_mg_day: float = 10.0,
    weeks: float = 8.0,
    age: int = 30,
    concurrent_antipsychotics: bool = False,
) -> ExposureProfile:
    """Build the structured exposure assessment for a planned course."""
    g = get_galactagogue(agent)
    cumulative = dose_mg_day * 7.0 * weeks

    d_tier = _duration_tier(weeks)
    a_tier = _age_tier(age)
    c_tier = _cumulative_tier(cumulative)
    tiers = {d_tier, a_tier, c_tier}
    if any(t.startswith("HIGH") for t in tiers):
        overall = "high"
    elif any(t.startswith("moderate") or t.startswith("elevated") for t in tiers):
        overall = "elevated"
    else:
        overall = "low"

    factors: List[str] = []
    if age >= 65:
        factors.append("advanced age (dominant TD risk factor)")
    if weeks > 12:
        factors.append("duration beyond the FDA 12-week recommendation")
    if cumulative > 6000:
        factors.append("high cumulative dose")
    if concurrent_antipsychotics:
        factors.append("concurrent dopamine antagonist (stacked risk)")
    if not factors:
        factors.append("no high-risk modifiers present in this profile")

    monitoring = [
        "Baseline: prolactin level + AIMS (Abnormal Involuntary Movement Scale) exam.",
        "Prolactin re-check ~week 2: confirm the dose actually reaches the target "
        "(this is the titration signal; the lowest dose that gets there is the goal).",
        "AIMS every 2-4 weeks during the course — this is a clinical exam; there "
        "is no blood test for TD. Any new involuntary movement = pause and reassess.",
        f"Hard stop-by cap: {g.name} labeling supports a 12-week ceiling; plan the "
        "taper, don't drift past it.",
        "AIMS again 1 and 3 months after stopping — TD can emerge or worsen "
        "post-exposure, so monitoring must outlive the drug.",
        "Once supply is established, step the dose DOWN to the minimum that "
        "maintains it (induction dose > maintenance dose).",
    ]

    notes = (
        "This is an exposure-framing tool, not a probability estimate: per-person "
        "TD incidence for short-course, low-dose, monitored use in young adults "
        "is not established. Population-average warnings are built on elderly, "
        "chronic-use cohorts and overstate the applicable risk for this profile — "
        "but the floor is not zero, rare short-course cases exist, and the failure "
        "mode (potentially irreversible movement disorder, no biomarker) is less "
        "forgiving than a lab-monitorable acute risk. The mitigation logic is "
        "identical to estradiol/VTE: lowest effective dose, shortest effective "
        "duration, structured monitoring — applied to a harsher failure mode."
    )

    return ExposureProfile(
        agent=g.name,
        dose_mg_day=dose_mg_day,
        weeks=weeks,
        age=age,
        cumulative_mg=cumulative,
        duration_tier=d_tier,
        age_tier=a_tier,
        cumulative_tier=c_tier,
        overall_tier=overall,
        factors=factors,
        monitoring=monitoring,
        notes=notes,
    )
