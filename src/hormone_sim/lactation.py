"""Lactation induction: protocol timeline + prolactin-driven milk production.

Models the **Newman–Goldfarb protocol** as used (off-label, case-report basis)
by trans women for induced lactation:

  **Phase 1 — Preparation.** Estradiol + anti-androgen, then progesterone
  added for the final weeks. This mimics pregnancy in two distinct steps:
    * **Ductal development** — driven by estradiol. Ducts are the scaffold.
    * **Lobuloalveolar (LA) development** — driven by progesterone, and it can
      only develop where the ductal scaffold already exists (P-gated by ducts).
      This is the "milk-making factory" the third trimester builds.

  The two are modeled as **separate states on separate timescales** because the
  published protocols differ in when P is added (≈1–3 months before the drop),
  and that timing is exactly what the model can interrogate.

  **Phase 2 — Transition ("the drop").** Abruptly stop estradiol and
  progesterone, simultaneously start the D2-antagonist galactagogue. This
  mimics delivery: losing progesterone removes the lactation brake and, with
  prolactin now high, triggers lactogenesis II. Tradeoff: dropping E2 also
  removes estrogenic support for breast tissue.

  **Phase 3 — Lactation.** Ongoing galactagogue + pumping/suckling. Milk
  supply ramps over weeks (lactogenesis maturation) toward a plateau set by
  glandular capacity.

Milk production model (phenomenological, population-typical):

    milk = Vmax * f_prolactin * capacity * f_suckling * maturation
    f_prolactin = P / (P + Km_P)
    capacity    = sqrt(ductal * la)        (both required; neither alone suffices)
    f_suckling  = sessions / (sessions + Km_sess)

``Vmax`` is the ceiling for a fully primed gland at saturating prolactin and
frequent stimulation (a fully lactating parent produces ~750–1000 mL/day).

**Epistemic status (read this):** there are NO controlled trials of lactation
induction in trans women. The landmark case report is Reisman & Goldstein,
*Induced lactation in a transgender woman*, N Engl J Med 2018 — one patient,
successful partial supply. This model is parameterized from postpartum cis
literature and medication-induced hyperprolactinemia data, then extrapolated.
``Vmax``, the Km values, the growth/maturation time constants, and the
involution rate are calibrated with LARGE uncertainty. Real case reports in
trans women describe partial supply (often tens-to-low-hundreds of mL/day),
not full lactation. Treat outputs as a *mechanistic reasoning tool*, not a
prediction.

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
    get_galactagogue,
    prolactin_response,
    PROLACTIN_BASELINE_NG_ML,
)
from .progestogens import get_progestogen


@dataclass(frozen=True)
class LactationModel:
    """Milk-production and glandular-development parameters.

    All calibrated with LARGE uncertainty — see module docstring. The two
    priming timescales are deliberately separate:

    * ``ductal_tau_weeks``   — E2-driven ductal scaffold (slow, months).
    * ``la_tau_weeks``       — progesterone-driven lobuloalveolar growth
      (faster, but gated by the ductal scaffold: LA cannot develop where
      ducts have not yet grown).
    """

    # Ceiling (mL/day) for a fully primed gland at saturating prolactin and
    # frequent stimulation. Full lactation ~750-1000 mL/day.
    milk_vmax_ml_day: float = 900.0

    # Prolactin (ng/mL) for half-maximal milk production. Lactation typically
    # needs prolactin well above baseline; term pregnancy ~200-400 ng/mL.
    milk_km_prolactin_ng_ml: float = 60.0

    # Ductal development: E2-driven, saturating. tau ~ months (pregnancy
    # ductal growth takes ~1st-2nd trimester; here 10 wk at reference E2).
    ductal_tau_weeks: float = 10.0

    # Lobuloalveolar development: progesterone-driven, gated by ducts.
    # Third-trimester-like P exposure of ~12 wk approaches LA ceiling.
    la_tau_weeks: float = 5.0

    # Estradiol effect on ductal growth is saturating with this EC50 (pg/mL).
    # 150 pg/mL (typical trans-feminine level) -> ~60% of maximal drive.
    ductal_e2_ec50_pgml: float = 100.0

    # Progesterone drive on lobuloalveolar growth: saturating in the
    # time-averaged P level (ng/mL) with this Km. Calibrated with large
    # uncertainty (see progestogens.py for why the mean level, not the peaks).
    la_p_km_ng_ml: float = 3.0

    # Lactogenesis II is triggered by progesterone WITHDRAWAL: maturation
    # progresses only when P has fallen below this fraction of its prep-phase
    # mean. This is why "the drop" matters — and why continuing P suppresses
    # supply even with a galactagogue on board.
    p_withdrawal_frac: float = 0.20

    # Estrogen BRAKE on established milk output. High-dose estrogen suppresses
    # lactation (historically used to do exactly that postpartum), so a
    # patient who keeps her HRT dose through the lactation phase pays a
    # supply cost. Brake fraction at saturating E2, and its EC50 (pg/mL).
    # Calibrated guess — flag accordingly.
    e2_brake_max: float = 0.60
    e2_brake_ec50_pgml: float = 200.0

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

    def f_prolactin(self, prolactin_ng_ml: float) -> float:
        p = max(0.0, float(prolactin_ng_ml))
        return p / (self.milk_km_prolactin_ng_ml + p)

    def f_suckling(self, sessions_per_day: float) -> float:
        s = max(0.0, float(sessions_per_day))
        return s / (self.suckling_km_sessions + s)

    def e2_drive(self, e2_pgml: float) -> float:
        """Saturating E2 drive for ductal growth, normalized to 1.0 at
        150 pg/mL (the reference trans-feminine level)."""
        e2 = max(0.0, float(e2_pgml))
        raw = e2 / (self.ductal_e2_ec50_pgml + e2)
        ref = 150.0 / (self.ductal_e2_ec50_pgml + 150.0)
        return raw / ref if ref > 0 else 0.0

    def involution_protection(self, sessions_per_day: float) -> float:
        """Fraction (0-1) of involution prevented by pumping frequency.

        Linear up to the recommended minimum, full protection at/above it —
        pumping maintains glandular structure; only *insufficient* pumping
        allows involution.
        """
        return min(1.0, max(0.0, float(sessions_per_day)) /
                   self.protection_sessions_per_day)

    def p_drive(self, p_level_ng_ml: float) -> float:
        """Saturating progesterone drive on LA development from the mean P
        level (ng/mL). The gland integrates exposure; see progestogens.py."""
        p = max(0.0, float(p_level_ng_ml))
        return p / (self.la_p_km_ng_ml + p)

    def p_withdrawn(self, p_level_ng_ml: float, prep_mean_ng_ml: float) -> bool:
        """Whether P has fallen far enough to permit lactogenesis II."""
        threshold = self.p_withdrawal_frac * max(prep_mean_ng_ml, 1e-9)
        return p_level_ng_ml <= threshold

    def e2_brake(self, e2_pgml: float) -> float:
        """Fraction (0-1) of milk output remaining given a maintained E2
        level during lactation (1.0 = no brake, at e2 = 0)."""
        e2 = max(0.0, float(e2_pgml))
        return 1.0 - self.e2_brake_max * (e2 / (self.e2_brake_ec50_pgml + e2))


@dataclass(frozen=True)
class LactationProtocol:
    """The Newman-Goldfarb-style protocol parameters.

    ``progesterone_start_weeks`` is when progesterone is added during the
    preparation phase (E2 + anti-androgen run from week 0 to ``prep_weeks``).
    Published protocols differ here (roughly the last 1-3 months); the model
    lets you interrogate that timing. ``None`` defaults to the last 12 weeks.
    Set it >= ``prep_weeks`` to run a P-free protocol (a genuinely open
    question — some protocols skip P).
    """

    prep_weeks: float = 30.0            # E2 + anti-androgen phase
    progesterone_start_weeks: Optional[float] = None  # None -> prep-12
    # Progesterone agent + dose: the P level now feeds a saturating LA drive
    # (PK-folded; see progestogens.py), so P dose-response is real in-model.
    progestogen: str = "micronized_progesterone"
    progesterone_dose_mg_day: float = 200.0
    # Transition shape. taper_weeks = 0 -> hard drop at prep_weeks (the
    # classic protocol). taper_weeks > 0 -> E2 (and P, if still on) decline
    # linearly to their post-drop levels over that many weeks AFTER the
    # drop; the galactagogue and pumping still start at prep_weeks. Tapering
    # delays the P-withdrawal trigger and thus "milk coming in".
    taper_weeks: float = 0.0
    # Estradiol level MAINTAINED during the lactation phase (pg/mL). 0 = full
    # drop (maximal supply, no tissue support). Keeping HRT-level E2 here
    # buys tissue support at a supply cost (the estrogen brake on lactation).
    post_drop_e2_pgml: float = 0.0
    galactagogue: str = "domperidone"
    galactagogue_dose_mg_day: float = 40.0
    pumping_sessions_per_day: float = 8.0
    # Estradiol level during preparation (pg/mL); drives ductal growth.
    prep_e2_pgml: float = 150.0

    def p_start(self) -> float:
        if self.progesterone_start_weeks is None:
            return max(0.0, self.prep_weeks - 12.0)
        return float(self.progesterone_start_weeks)

    def summary(self) -> str:
        p_on = self.p_start()
        p_txt = f"{p_on:g}-{self.prep_weeks:g} wk" if p_on < self.prep_weeks else "none"
        taper_txt = (f"{self.taper_weeks:g} wk taper"
                     if self.taper_weeks > 0 else "hard drop")
        e2_txt = (f"{self.post_drop_e2_pgml:g} pg/mL maintained"
                  if self.post_drop_e2_pgml > 0 else "full drop")
        return (
            f"Protocol: Newman-Goldfarb-style induction\n"
            f"  Prep: {self.prep_weeks:g} wk on E2 ({self.prep_e2_pgml:g} pg/mL) "
            f"+ anti-androgen\n"
            f"  Progesterone: {p_txt} "
            f"({self.progesterone_dose_mg_day:g} mg/d "
            f"{self.progestogen.replace('_', ' ')})\n"
            f"  Transition: {taper_txt}; post-drop E2: {e2_txt}\n"
            f"  Galactagogue: {self.galactagogue} "
            f"{self.galactagogue_dose_mg_day:g} mg/day\n"
            f"  Pumping: {self.pumping_sessions_per_day:g} sessions/day"
        )


@dataclass
class LactationResult:
    """Outcome of a simulated lactation induction."""

    protocol: LactationProtocol
    weeks: np.ndarray
    ductal: np.ndarray            # 0-1 E2-driven ductal scaffold
    lobuloalveolar: np.ndarray    # 0-1 P-driven milk-secretory units
    priming: np.ndarray           # combined capacity index = sqrt(ductal*LA)
    maturation: np.ndarray        # 0-1 lactogenesis maturation (milk coming in)
    prolactin: np.ndarray         # ng/mL (weekly average)
    p_level_ng_ml: np.ndarray     # mean progesterone level (ng/mL)
    e2_level_pgml: np.ndarray     # estradiol level (pg/mL)
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
        drop_i = int(np.argmax(self.weeks >= self.protocol.prep_weeks))
        lines = [
            self.protocol.summary(),
            f"  Prolactin on galactagogue (steady state): "
            f"{prolactin_response(g, self.protocol.galactagogue_dose_mg_day):.0f} ng/mL",
            f"  At the drop: ductal {self.ductal[drop_i]:.2f}, "
            f"lobuloalveolar {self.lobuloalveolar[drop_i]:.2f} "
            f"(capacity {self.priming[drop_i]:.2f})",
            f"  Milk: peak {self.peak_milk():.0f} mL/day, "
            f"steady-state (last 8 wk) {ss:.0f} mL/day",
            f"  Meaningful supply (>50 mL/day): "
            f"{f'{first:.1f} weeks after the drop' if first is not None else 'not reached'}",
            "",
            f"  Model: milk = Vmax {m.milk_vmax_ml_day:.0f} * f(P; Km "
            f"{m.milk_km_prolactin_ng_ml:.0f}) * sqrt(ductal*LA) * f(sessions; "
            f"Km {m.suckling_km_sessions:.0f}) * maturation(tau "
            f"{m.maturation_tau_weeks:.0f} wk)",
            "  WARNING: no controlled trials in trans women; calibrated from "
            "postpartum cis literature (case reports report partial supply). "
            "Not medical advice.",
        ]
        return "\n".join(lines)


def simulate_lactation(
    protocol: LactationProtocol,
    total_weeks: Optional[float] = None,
    dt_weeks: float = 0.25,
    model: Optional[LactationModel] = None,
    individual_prolactin_ng_ml: Optional[float] = None,
) -> LactationResult:
    """Simulate the full protocol timeline.

    Prep: ductal grows with E2; lobuloalveolar grows with the *mean*
    progesterone level (PK-folded via progestogens.py), gated by the ductal
    scaffold. Transition: hard drop by default, or a linear taper of E2 (and
    P) over ``taper_weeks``. Lactation: galactagogue raises prolactin, but
    maturation (lactogenesis II) waits for progesterone WITHDRAWAL; pumping
    maintains (never rebuilds) glandular structure; any E2 maintained through
    this phase applies the estrogen brake to output.

    ``individual_prolactin_ng_ml`` overrides the agent's population dose-
    response (used by the personalization/sensitivity layers to inject a
    sampled individual's prolactin level).
    """
    model = model or LactationModel()
    g = get_galactagogue(protocol.galactagogue)
    prog = get_progestogen(protocol.progestogen)

    if total_weeks is None:
        # Prep + a year of lactation phase by default.
        total_weeks = protocol.prep_weeks + 52.0

    n = int(round(total_weeks / dt_weeks)) + 1
    weeks = np.arange(n) * dt_weeks

    ductal = np.zeros(n)
    la = np.zeros(n)
    capacity = np.zeros(n)
    maturation = np.zeros(n)
    prolactin = np.zeros(n)
    p_level = np.zeros(n)
    e2_level = np.zeros(n)
    milk = np.zeros(n)
    phase = np.zeros(n, dtype=int)

    # Ductal / LA growth constants (per week, at reference drive).
    k_d = 1.0 / model.ductal_tau_weeks
    k_la = 1.0 / model.la_tau_weeks

    # --- hormone level trajectories (mean levels on the weekly grid) --------
    p_start = protocol.p_start()
    p_end = protocol.prep_weeks
    p_on_mean = prog.mean_conc_ng_ml(protocol.progesterone_dose_mg_day)
    taper = max(0.0, float(protocol.taper_weeks))
    taper_end = protocol.prep_weeks + taper

    for i, w in enumerate(weeks):
        in_lactation = w >= protocol.prep_weeks
        phase[i] = 1 if in_lactation else 0

        # Progesterone: on during the P window; linear taper after the drop
        # if tapering; zero otherwise.
        if w < p_start or w >= p_end + taper:
            p_level[i] = 0.0
        elif w < p_end:
            p_level[i] = p_on_mean
        else:  # taper window: linear decline prep_mean -> 0
            frac = (w - p_end) / taper if taper > 0 else 1.0
            p_level[i] = p_on_mean * max(0.0, 1.0 - frac)

        # Estradiol: prep level until the drop; linear taper to the
        # post-drop level if tapering; maintenance level afterwards.
        if w < protocol.prep_weeks:
            e2_level[i] = protocol.prep_e2_pgml
        elif w < taper_end and taper > 0:
            frac = (w - protocol.prep_weeks) / taper
            e2_level[i] = (protocol.prep_e2_pgml +
                           (protocol.post_drop_e2_pgml - protocol.prep_e2_pgml) * frac)
        else:
            e2_level[i] = protocol.post_drop_e2_pgml

    # Prolactin: population dose-response, or a sampled individual's level.
    if individual_prolactin_ng_ml is not None:
        p_drug = float(individual_prolactin_ng_ml)
    else:
        p_drug = prolactin_response(g, protocol.galactagogue_dose_mg_day)

    # Maturation gates.
    p_threshold = model.maturation_prolactin_threshold_ng_ml
    tau_m = model.maturation_tau_weeks

    for i, w in enumerate(weeks):
        in_lactation = w >= protocol.prep_weeks

        if not in_lactation:
            # Prep: ductal grows with the current E2 level; LA grows with
            # the mean P drive, gated by the ductal scaffold.
            d_prev = ductal[i - 1] if i > 0 else 0.0
            la_prev = la[i - 1] if i > 0 else 0.0
            ductal[i] = min(1.0, d_prev + k_d * model.e2_drive(e2_level[i])
                            * (1.0 - d_prev) * dt_weeks)
            p_on = 1.0 if p_level[i] > 0 else 0.0
            drive = model.p_drive(p_level[i]) * p_on
            la[i] = min(1.0, la_prev + k_la * drive * ductal[i]
                        * (1.0 - la_prev) * dt_weeks)
            maturation[i] = 0.0
            prolactin[i] = PROLACTIN_BASELINE_NG_ML
            milk[i] = 0.0
        else:
            # Lactation: pumping MAINTAINS glandular structure (prevents
            # involution) but does not rebuild it — ductal/LA architecture
            # comes from E2+progesterone during the prep phase.
            d_prev = ductal[i - 1] if i > 0 else 1.0
            la_prev = la[i - 1] if i > 0 else 1.0
            protection = model.involution_protection(protocol.pumping_sessions_per_day)
            decay = model.involution_rate_per_week * (1.0 - protection) * dt_weeks
            ductal[i] = max(0.0, d_prev - decay)
            la[i] = max(0.0, la_prev - decay)

            # Prolactin rises (galactagogue); maturation follows with a lag,
            # but ONLY once progesterone has withdrawn (lactogenesis II).
            prolactin[i] = p_drug
            prev_m = maturation[i - 1] if i > 0 else 0.0
            gates_open = (prolactin[i] > p_threshold and
                          model.p_withdrawn(p_level[i], p_on_mean))
            if gates_open:
                maturation[i] = min(1.0, prev_m + (1.0 - prev_m) * dt_weeks / tau_m)
            else:
                maturation[i] = prev_m

            f_p = model.f_prolactin(prolactin[i])
            f_s = model.f_suckling(protocol.pumping_sessions_per_day)
            brake = model.e2_brake(e2_level[i])
            milk[i] = (model.milk_vmax_ml_day * f_p * f_s * brake * maturation[i])

        capacity[i] = float(np.sqrt(max(0.0, ductal[i] * la[i])))

    milk *= capacity

    return LactationResult(
        protocol=protocol,
        weeks=weeks,
        ductal=ductal,
        lobuloalveolar=la,
        priming=capacity,
        maturation=maturation,
        prolactin=prolactin,
        p_level_ng_ml=p_level,
        e2_level_pgml=e2_level,
        milk_ml_day=milk,
        phase=phase,
    )
