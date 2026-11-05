"""Lactation induction: protocol timeline + prolactin-driven milk production.

Models the **Newman–Goldfarb protocol** as used (off-label, case-report basis)
by trans women for induced lactation:

  **Phase 1 — Preparation.** Estradiol + anti-androgen + progesterone for
  ~6–9 months. This mimics pregnancy: E2 drives ductal growth, progesterone
  drives lobuloalveolar development, and together they raise prolactin
  modestly. The output is a **priming factor** — how far the glandular
  tissue has been brought toward a milk-capable state.

  **Phase 2 — Transition ("the drop").** Abruptly stop estradiol and
  progestogen, simultaneously start the D2-antagonist galactagogue. This
  mimics delivery: losing placental progesterone removes the lactation brake
  and, with prolactin now high, triggers lactogenesis II. Note the tradeoff:
  dropping E2 also removes the estrogenic support for breast tissue.

  **Phase 3 — Lactation.** Ongoing galactagogue + pumping/suckling.
  Milk supply ramps over days-to-weeks as the gland responds to prolactin
  and frequent stimulation.

Milk production model (phenomenological, population-typical):

    milk = Vmax * f_prolactin * priming * f_suckling
    f_prolactin = P / (P + Km_P)           (saturating, Km ~ lactation threshold)
    f_suckling  = sessions / (sessions + Km_sess)

``Vmax`` is the ceiling for a well-primed gland at saturating prolactin and
frequent stimulation (a fully lactating parent produces ~750–1000 mL/day).

**Epistemic status (read this):** there are NO controlled trials of lactation
induction in trans women. This model is parameterized from postpartum cis
literature and medication-induced hyperprolactinemia data, then extrapolated.
``Vmax``, ``Km_P``, priming rate, and involution rate are calibrated with LARGE
uncertainty. Real case reports in trans women describe partial supply
(often tens-to-low-hundreds of mL/day), not full lactation. Treat outputs as
a *mechanistic reasoning tool*, not a prediction.

Not medical advice: domperidone is not FDA-approved (QT risk), metoclopramide
carries a boxed warning (tardive dyskinesia). Lactation induction should be
run with a clinician monitoring prolactin, cardiac status, and
neuropsychiatric effects.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .galactagogues import (
    Galactagogue,
    get_galactagogue,
    prolactin_response,
    PROLACTIN_BASELINE_NG_ML,
    PROLACTIN_LACTATION_REFERENCE_NG_ML,
)


@dataclass(frozen=True)
class LactationModel:
    """Milk-production parameters (calibrated, LARGE uncertainty)."""

    # Ceiling (mL/day) for a fully primed gland at saturating prolactin and
    # frequent stimulation. Full lactation ~750-1000 mL/day.
    milk_vmax_ml_day: float = 900.0

    # Prolactin (ng/mL) for half-maximal milk production. Lactation typically
    # needs prolactin well above baseline; term pregnancy ~200-400 ng/mL.
    milk_km_prolactin_ng_ml: float = 60.0

    # Priming: saturating growth during preparation (E2+P phase). tau in weeks —
    # ~6-9 months of prep gets priming to ~0.9-0.95.
    priming_tau_weeks: float = 10.0

    # Lactogenesis maturation: the gland takes time to respond to elevated
    # prolactin ("milk coming in"). tau_maturation sets that ramp (weeks).
    # Cis postpartum: 2-5 days; trans-feminine case reports: ~2-6 weeks.
    maturation_tau_weeks: float = 3.0

    # Prolactin threshold (ng/mL) above which maturation progresses. Below it
    # the gland does not further mature.
    maturation_prolactin_threshold_ng_ml: float = 40.0

    # Suckling/pumping: sessions per day for half-maximal supply signal.
    # Newman-Goldfarb recommends 8-12 sessions/day.
    suckling_km_sessions: float = 8.0

    # Stimulation sufficient to fully prevent involution (sessions/day).
    # Clinical guidance: 8-12 sessions/day maintains supply; below ~6,
    # involution/supply-loss sets in.
    protection_sessions_per_day: float = 6.0

    # Priming decay (involution) if not stimulated after the drop, per week.
    # Weaning-level involution happens over weeks-to-months.
    involution_rate_per_week: float = 0.05

    def involution_protection(self, sessions_per_day: float) -> float:
        """Fraction (0-1) of involution prevented by pumping frequency.

        Linear up to the recommended minimum, full protection at/above it —
        pumping maintains glandular structure; only *insufficient* pumping
        allows involution.
        """
        return min(1.0, max(0.0, float(sessions_per_day)) /
                   self.protection_sessions_per_day)

    def f_prolactin(self, prolactin_ng_ml: float) -> float:
        p = max(0.0, float(prolactin_ng_ml))
        return p / (self.milk_km_prolactin_ng_ml + p)

    def f_suckling(self, sessions_per_day: float) -> float:
        s = max(0.0, float(sessions_per_day))
        return s / (self.suckling_km_sessions + s)


@dataclass(frozen=True)
class LactationProtocol:
    """The Newman-Goldfarb-style protocol parameters."""

    prep_weeks: float = 30.0            # ~7 months of E2+P priming
    galactagogue: str = "domperidone"
    galactagogue_dose_mg_day: float = 40.0
    pumping_sessions_per_day: float = 8.0
    # Estradiol level during preparation (pg/mL), for context/priming only.
    prep_e2_pgml: float = 150.0

    def summary(self) -> str:
        return (
            f"Protocol: Newman-Goldfarb-style induction\n"
            f"  Prep: {self.prep_weeks:g} weeks on E2 ({self.prep_e2_pgml:g} pg/mL) "
            f"+ anti-androgen + progesterone\n"
            f"  Galactagogue: {self.galactagogue} {self.galactagogue_dose_mg_day:g} mg/day\n"
            f"  Pumping: {self.pumping_sessions_per_day:g} sessions/day"
        )


@dataclass
class LactationResult:
    """Outcome of a simulated lactation induction."""

    protocol: LactationProtocol
    weeks: np.ndarray
    priming: np.ndarray           # 0-1 glandular priming index
    maturation: np.ndarray        # 0-1 lactogenesis maturation (milk coming in)
    prolactin: np.ndarray         # ng/mL (weekly average)
    milk_ml_day: np.ndarray       # mL/day
    phase: np.ndarray             # 0=prep, 1=lactation

    def steady_state_milk(self, tail_weeks: float = 8.0) -> float:
        """Mean milk (mL/day) over the last ``tail_weeks``."""
        dt = self.weeks[1] - self.weeks[0] if self.weeks.size > 1 else 1.0
        n_tail = max(1, int(round(tail_weeks / dt)))
        return float(self.milk_ml_day[-n_tail:].mean())

    def peak_milk(self) -> float:
        return float(self.milk_ml_day.max())

    def weeks_to_first_milk(self, threshold_ml_day: float = 50.0) -> Optional[float]:
        """Weeks after the drop when milk first exceeds ``threshold_ml_day``
        (a 'meaningful supply' threshold, not the very first drop)."""
        lactation_start = self.protocol.prep_weeks
        for i, w in enumerate(self.weeks):
            if w >= lactation_start and self.milk_ml_day[i] >= threshold_ml_day:
                return float(w - lactation_start)
        return None

    def report(self) -> str:
        m = LactationModel()
        ss = self.steady_state_milk()
        first = self.weeks_to_first_milk()
        g = get_galactagogue(self.protocol.galactagogue)
        lines = [
            self.protocol.summary(),
            f"  Prolactin on galactagogue (steady state): "
            f"{prolactin_response(g, self.protocol.galactagogue_dose_mg_day):.0f} ng/mL",
            f"  Priming at drop: "
            f"{self.priming[np.argmax(self.weeks >= self.protocol.prep_weeks)]:.2f}",
            f"  Milk: peak {self.peak_milk():.0f} mL/day, "
            f"steady-state (last 8 wk) {ss:.0f} mL/day",
            f"  Meaningful supply (>50 mL/day): "
            f"{f'{first:.1f} weeks after the drop' if first is not None else 'not reached'}",
            "",
            f"  Model: milk = Vmax {m.milk_vmax_ml_day:.0f} * f(P; Km "
            f"{m.milk_km_prolactin_ng_ml:.0f}) * priming * f(sessions; Km "
            f"{m.suckling_km_sessions:.0f}) * maturation(tau "
            f"{m.maturation_tau_weeks:.0f} wk)",
            "  WARNING: no controlled trials in trans women; calibrated from "
            "postpartum cis literature. Partial supply is the reported norm, "
            "not full lactation. Not medical advice.",
        ]
        return "\n".join(lines)


def simulate_lactation(
    protocol: LactationProtocol,
    total_weeks: Optional[float] = None,
    dt_weeks: float = 0.25,
    model: Optional[LactationModel] = None,
) -> LactationResult:
    """Simulate the full protocol timeline.

    Prep phase: priming grows saturatingly (no galactagogue).
    Lactation phase: galactagogue at steady state + pumping; priming is
    maintained by stimulation (decays only if sessions ~ 0).
    """
    model = model or LactationModel()
    g = get_galactagogue(protocol.galactagogue)

    if total_weeks is None:
        # Prep + a year of lactation phase by default.
        total_weeks = protocol.prep_weeks + 52.0

    n = int(round(total_weeks / dt_weeks)) + 1
    weeks = np.arange(n) * dt_weeks

    priming = np.zeros(n)
    maturation = np.zeros(n)
    prolactin = np.zeros(n)
    milk = np.zeros(n)
    phase = np.zeros(n, dtype=int)

    # Prep-phase priming: saturating approach to 1.0 with tau.
    tau = model.priming_tau_weeks

    # Prolactin on the galactagogue (steady-state, dose-based) — the PK ramp
    # is fast (days) relative to the weekly grid, so use the steady state.
    p_drug = prolactin_response(g, protocol.galactagogue_dose_mg_day)

    # Maturation grows only while prolactin exceeds the threshold.
    p_threshold = model.maturation_prolactin_threshold_ng_ml
    tau_m = model.maturation_tau_weeks

    for i, w in enumerate(weeks):
        in_lactation = w >= protocol.prep_weeks
        phase[i] = 1 if in_lactation else 0

        if not in_lactation:
            # Prep: priming grows; prolactin stays near baseline.
            priming[i] = 1.0 - np.exp(-w / tau)
            prolactin[i] = PROLACTIN_BASELINE_NG_ML
            maturation[i] = 0.0
            milk[i] = 0.0
        else:
            # Lactation: pumping MAINTAINS the gland's development (prevents
            # involution) but does not rebuild it — ductal/lobuloalveolar
            # architecture comes from E2+progesterone during the prep phase.
            # So priming can only decay here, more slowly when stimulated.
            prev = priming[i - 1] if i > 0 else 1.0
            stim_factor = model.f_suckling(protocol.pumping_sessions_per_day)
            protection = model.involution_protection(protocol.pumping_sessions_per_day)
            decay = model.involution_rate_per_week * (1.0 - protection)
            priming[i] = max(0.0, prev - decay * dt_weeks)

            # Prolactin rises; maturation follows with a lag when above threshold.
            prolactin[i] = p_drug
            prev_m = maturation[i - 1] if i > 0 else 0.0
            if prolactin[i] > p_threshold:
                maturation[i] = min(1.0, prev_m + (1.0 - prev_m) * dt_weeks / tau_m)
            else:
                maturation[i] = prev_m

            f_p = model.f_prolactin(prolactin[i])
            f_s = stim_factor
            milk[i] = model.milk_vmax_ml_day * f_p * priming[i] * f_s * maturation[i]

    return LactationResult(
        protocol=protocol,
        weeks=weeks,
        priming=priming,
        maturation=maturation,
        prolactin=prolactin,
        milk_ml_day=milk,
        phase=phase,
    )
