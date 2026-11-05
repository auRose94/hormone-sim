"""Progestogen pharmacokinetics for the lactation-induction prep phase.

Folds real progesterone PK into the lactation model, replacing the crude
on/off switch the P phase previously used.

**Micronized oral progesterone has wild peaks** — this is the pharmacological
reality the user asked to see modeled:

* t½ ~6 h with once-daily (bedtime) dosing → serum oscillates between a
  near-zero trough and a large peak every day.
* Low oral bioavailability (~10%, first-pass) but a very large volume of
  distribution (~10 L/kg) — the combination that produces mean levels around
  ~10 ng/mL at 200 mg/day while Cmax reaches ~30 ng/mL.
* Consequence for the model: lobuloalveolar development **integrates**
  exposure (the gland responds to sustained P, not to individual peaks), so
  the LA drive uses the *time-averaged* level. The peaks are not irrelevant —
  they drive the sedating side effects and they make the **withdrawal**
  abrupt — but modeling intra-day oscillation on a weekly simulation grid
  would be fake precision. The PK module therefore provides mean level,
  expected peak, and peak/mean ratio; the lactation model uses the mean
  trajectory and the documented peaks appear in reports.

**Withdrawal is the other half of the physiology.** Progesterone actively
*brakes* milk secretion; lactogenesis II is triggered by P **withdrawal**
(prolactin being high). The lactation model gates maturation on P falling
below a fraction of its prep level — which is why "the drop" matters and why
continuing P into the lactation phase suppresses supply even with a
galactagogue on board.

Parameter provenance: micronized progesterone PK is well characterized
(peak 1-3 h post-dose, t½ ~6 h, F ~8-15%, Vd large) — see
docs/references.md. The *mapping* from P exposure to LA development (Km of
the drive curve) is calibrated with large uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np


@dataclass(frozen=True)
class Progestogen:
    """A progestogen's PK for the prep-phase P exposure."""

    name: str
    half_life_h: float
    bioavailability: float
    vd_l_kg: float
    typical_dose_mg_day: float
    citation: str
    note: str = ""

    def body_weight_kg(self) -> float:
        return 70.0

    def clearance_l_per_h(self) -> float:
        """Cl = ln2 * Vd / t1/2."""
        return (0.6931471805599453 * self.vd_l_kg * self.body_weight_kg()
                / self.half_life_h)

    def mean_conc_ng_ml(self, dose_mg_day: float) -> float:
        """Time-averaged serum level (ng/mL) at a given daily oral dose.

        ng/mL == µg/L: mean input (µg/h) / clearance (L/h).
        """
        input_ug_per_h = dose_mg_day * 1e3 * self.bioavailability / 24.0
        return input_ug_per_h / self.clearance_l_per_h()

    def peak_conc_ng_ml(self, dose_mg_day: float,
                        dosing_interval_h: float = 24.0) -> float:
        """Peak serum level (ng/mL) at steady state with q24h dosing.

        Bolus approximation with accumulation factor
        1 / (1 - exp(-k * tau)). Real micronized-P peaks are somewhat lower
        than a bolus predicts (absorption takes 1-3 h) — treat as an upper
        estimate.
        """
        k = 0.6931471805599453 / self.half_life_h
        acc = 1.0 / (1.0 - np.exp(-k * dosing_interval_h))
        vd_ml = self.vd_l_kg * self.body_weight_kg() * 1e3
        single = (dose_mg_day * 1e6 * self.bioavailability) / vd_ml  # ng/mL
        return single * acc

    def peak_to_mean_ratio(self, dose_mg_day: float) -> float:
        """How 'wild' the peaks are at this dose (Cmax / Cmean)."""
        mean = self.mean_conc_ng_ml(dose_mg_day)
        if mean <= 0:
            return 0.0
        return self.peak_conc_ng_ml(dose_mg_day) / mean


# --- Registry ------------------------------------------------------------------
# Micronized progesterone is the agent used in trans-feminine induction
# protocols. Only well-characterized agents are registered; a custom
# constructor exists for exploring alternatives without fabricating provenance.

MICRONIZED_PROGESTERONE = Progestogen(
    name="micronized_progesterone",
    half_life_h=6.0,
    bioavailability=0.10,
    vd_l_kg=10.0,
    typical_dose_mg_day=200.0,
    citation="Micronized progesterone PK: peak 1-3 h, t1/2 ~6 h, F ~8-15%, "
    "Vd ~10 L/kg (StatPearls 'Progesterone'; Simon et al. 1993 dose-ranging). "
    "200 mg/day -> mean ~10 ng/mL, Cmax ~25-30 ng/mL.",
    note="The 'wild peaks' agent: q24h dosing with t1/2 6 h gives a ~3x "
    "peak/mean swing and a near-zero trough. Bedtime dosing is standard "
    "because of sedation. LA development integrates exposure; the peaks "
    "matter for side effects and for how abruptly withdrawal lands.",
)

PROGESTOGENS: Dict[str, Progestogen] = {
    "micronized_progesterone": MICRONIZED_PROGESTERONE,
}


def get_progestogen(name: str) -> Progestogen:
    if name not in PROGESTOGENS:
        raise KeyError(
            f"Unknown progestogen {name!r}. Known: {sorted(PROGESTOGENS)}. "
            "Use custom_progestogen() to explore an alternative explicitly."
        )
    return PROGESTOGENS[name]


def custom_progestogen(
    name: str,
    half_life_h: float,
    bioavailability: float,
    vd_l_kg: float,
    typical_dose_mg_day: float,
    note: str = "",
) -> Progestogen:
    """Construct a progestogen with user-specified PK (for exploring
    alternatives — e.g. a long-half-life synthetic — without pretending the
    registry characterized it)."""
    return Progestogen(
        name=name, half_life_h=half_life_h, bioavailability=bioavailability,
        vd_l_kg=vd_l_kg, typical_dose_mg_day=typical_dose_mg_day,
        citation="user-specified PK (custom progestogen)",
        note=note or "Custom progestogen; PK supplied by caller, not literature.",
    )
