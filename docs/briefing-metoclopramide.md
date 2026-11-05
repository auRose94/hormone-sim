# Briefing: low-dose, short-course metoclopramide for induced lactation in trans women

**Prepared for:** clinician discussion (trans-inclusive primary care / endocrinology)
**Status of this document:** a mechanistic modeling analysis + literature synthesis.
It is **not** a recommendation, and the underlying model is **not** validated — see
"Evidence quality" before anything else.

---

## 1. Evidence quality — read this first

* There are **no controlled trials** of pharmacologic lactation induction in
  trans women. The closest datapoint to a trans-feminine validation is **one
  case report**: Reisman T, Goldstein Z. *Induced lactation in a transgender
  woman.* N Engl J Med 2018;378 — successful **partial** supply.
* Galactagogue PK and prolactin dose-response come from **postpartum cis women**
  and medication-induced hyperprolactinemia literature, extrapolated here to
  trans-feminine physiology (on estradiol, no pregnancy priming).
* The model (`hormone-sim`, open source, reproducible; see repo) is
  population-typical with **large, explicitly flagged uncertainty**. Its outputs
  are order-of-magnitude reasoning tools, and the sensitivity bands are the
  honest product — not the point estimates.

## 2. The question this analysis addresses

Whether a **low-dose, short-course** metoclopramide protocol can occupy a
useful/danger-separated region: enough prolactin elevation for lactogenesis,
with TD (tardive dyskinesia) exposure held near the floor. The patient
hypothesis motivating this: *efficacy saturates with dose; TD exposure does
not — it accumulates linearly with dose × time.* If true, a "knee" exists.

## 3. What the model says (n=150 paired samples per dose; identical protocol:
30 wk E2/anti-androgen prep with P in final 12 wk, drop, galactagogue + 8 pumping sessions/day)

| Dose | Prolactin (ng/mL) | Milk p10 / p50 / p90 (mL/day) | % of modeled bodies ≥100 mL/day |
|---|---|---|---|
| 5 mg/d | 37 | 0 / 0 / 195 | 40% |
| 10 mg/d | 55 | 0 / 173 / 268 | 81% |
| 15 mg/d | 68 | 138 / 204 / 302 | 94% |
| 20 mg/d | 79 | 164 / 222 / 320 | 98% |
| 30 mg/d | 94 | 179 / 239 / 341 | 99% |
| 45 mg/d | 108 | 189 / 254 / 358 | 99% |

Findings, stated carefully:

1. **The knee is at ~15–20 mg/day.** Past 20 mg/d, each added mg buys ~1–2
   mL/day of median supply; below 10 mg/d, a meaningful fraction of modeled
   physiologies fail to respond at all (p10 = 0 at 5–10 mg/d).
2. **The low-dose zone is heterogeneous, not uniformly weak.** At 10 mg/d,
   ~81% of modeled bodies reach 100 mL/day but the rest do not — driven by
   (unknown, sampled) spread in prolactin sensitivity and the gland's
   maturation threshold. **This is an argument for titration-to-effect, not
   fixed low dosing:** find the individual's minimum effective dose.
3. **Exposure asymmetry:** at 15 mg/d × 12 wk the cumulative exposure is
   1,260 mg; the chronic-use profiles behind population TD warnings
   (e.g., 40 mg/d × 30 wk, elderly) are ~8,400 mg. Same drug, ~7× exposure —
   the applicable-risk argument for young short-course patients is
   pharmacologically legitimate, with the caveat in §5.

## 4. Proposed protocol shape (if a study were run)

* **Baseline:** prolactin, AIMS (Abnormal Involuntary Movement Scale), pregnancy/lactation history, concurrent medications.
* **Induction:** metoclopramide 10 mg/day (divided) → prolactin check at wk 2 →
  titrate in 5 mg/day steps to the **lowest dose reaching target** (e.g.,
  ~60–80 ng/mL, the model's half-maximal-milk zone), ceiling 20 mg/day.
* **Duration:** hard 12-week ceiling per FDA labeling; pre-planned taper.
* **Maintenance:** once supply established, step down to minimum maintaining dose.
* **Monitoring:** AIMS every 2–4 weeks during, and at 1 and 3 months
  **post-stop** (TD can emerge after exposure ends). No biomarker exists;
  the exam is the monitor.
* **Endpoints:** minimal effective dose distribution; supply trajectory vs
  model bands (pre-registered as expected values); TD incidence (expected 0
  at feasible n — report as such); prolactin achieved vs dose.

This is deliberately an **n-of-1 titration core wrapped in a case series** —
the cheapest design that produces individually meaningful data while building
toward the cohort that doesn't exist.

## 5. What the model cannot say (and the honest risk statement)

* It cannot give a per-person TD probability — that number does not exist for
  this exposure profile. The risk tool (`hormone-sim risk`) outputs *exposure
  framing* (tiers from duration/age/cumulative dose), explicitly not prophecy.
* The asymmetry vs thrombotic risk stands: VTE is acute, lab-monitorable,
  reversible; TD is cumulative, exam-monitorable only, potentially
  irreversible, and can emerge post-exposure. Population warnings **overstate**
  the applicable risk for a young short-course patient and the floor is
  **not zero**. Both facts, simultaneously.
* Prolactin dose-response in trans-feminine physiology is itself unmeasured
  (the model borrows it from cis populations); per-patient prolactin
  measurement is therefore not optional — it is the primary titration signal.

## 6. Pointers

Model, parameters, and provenance: `hormone-sim` repository —
`docs/references.md` (all parameter citations, flagged calibrated vs measured),
`docs/governance.md` (scope/ethics), `src/hormone_sim/risk_exposure.py`
(risk framing), `src/hormone_sim/sensitivity.py` (bands and sweep).
Key external sources: Reisman & Goldstein NEJM 2018 (n=1 case);
Endocrine Society CPG 2017; WPATH SOC-8 (2022); FDA metoclopramide labeling
(boxed warning, 12-week guidance); galactagogue reviews (Sankar 2017;
OSHERC domperidone 2018).
