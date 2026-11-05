# Roadmap: trans-masculine / masculinization wave (hormone-sim v0.8+)

**Status: planned, not built.** This document exists so the next wave starts
from a plan instead of a blank file — and so anyone picking the project up
later (including other trans researchers, which is the point of sharing it)
can see the intended architecture before writing code.

The trans-feminine work established a reusable architecture: registry-style
compound PK, Hill-function PD, population bands, heterogeneity sweeps,
risk-exposure framing, and the personalization (importance-sampling) layer.
The masculinization wave is largely **new registries + new target ranges +
new PD outcomes** on the same chassis — plus honest acknowledgment that the
trans-masculine evidence base has its own gaps, some worse than the
trans-feminine ones.

## 1. Testosterone PK (new registry: `androgens.py`)

| Agent | Route | PK notes to model |
|---|---|---|
| Testosterone cypionate / enanthate | IM / SubQ | Depot ester, t½ ~8–11 d (cypionate), classic 100–200 mg/1–2 wk; large peak–trough swings — the trans-masc analog of "wild peaks," and the trough is where dysphoric symptoms often land |
| Testosterone enanthate (shorter interval) | SubQ | Weekly/biweekly small volumes; smoother curves at the same weekly dose |
| Testosterone gel / cream | transdermal | Daily zero-order-ish; steady levels; lower peak; absorption varies by site/skin |
| Testosterone undecanoate | oral / IM | Oral: with-fat requirement, variable F; IM: very long depot (10–14 wk intervals) |

Same design rule as estradiol: clearance/Vd are molecule properties; route
differences live in the input function + calibrated delivery factors. The
galactagogue PK pattern (one-compartment, repeated dosing, accumulation) is
directly reusable for the esters.

## 2. Testosterone PD / target ranges (extend `params.py` + new `pd` branch)

* **Target:** male-typical T, 300–1000 ng/dL trough-adjusted (Endocrine
  Society CPG); the model should track *trough* vs *peak* distinctly since
  trough drives symptom patterns.
* **Suppression direction reversed:** T (via aromatization to estradiol and
  direct androgen feedback) suppresses LH/FSH → ovarian E2/progesterone
  production falls. The existing Hill-feedback machinery works with new
  anchors: baseline AMAB-female E2 (~30–80 pg/mL pre-GAHT... cite properly
  at build time), suppression EC50 to calibrate to "cycles cease by ~3–6 mo
  on stable T."
* **5α-DHT** (skin/hair/voice-relevant metabolite) as a second PD state if
  justifiable from data — or explicitly omitted with a note, rather than
  hand-waved.

## 3. Masculinization PD outcomes (new `masculinization.py`, mirrors `breast_development.py`)

Cumulative-exposure-driven, saturating, population-typical — same honest
"genetics dominate" framing as the breast model:

* **Voice drop** (pitch decline over ~6–24 mo; partial, not always to
  cis-male range — voice training is a separate modality, link out)
* **Body/facial hair** (Tanner-like staging over years)
* **Fat redistribution** + muscle mass (T-anchored, slower)
* **Menstrual cessation** (usually 2–6 mo — the one reasonably well-tracked
  outcome, make it the calibration anchor)
* **Bottom growth** (clitoromegaly, early and T-dose-sensitive)

## 4. Risk framing for trans-masc (the analog of the VTE/TD discussion)

* **Erythrocytosis** (Hct > 50–52%): the *dose-and-route-dependent,
  lab-monitorable* risk — this is the trans-masc case where the user's
  "monitorable vs not" asymmetry argument lands the *other* way: it IS a
  blood-monitorable risk with dose/route heterogeneity. Risk-exposure module
  analogous to `risk_exposure.py`, with the honest note that polycythemia
  management (dose reduction, donation) is well described.
* **Lipids / cardiovascular** — weaker evidence; frame as monitoring, not
  prophecy.
* **Atrophic vaginitis / UTI** from E2 suppression — common, manageable,
  rarely modeled anywhere.
* **Mental-health trough effects** — report as a known patient-reported
  pattern tied to peak–trough PK shape; the model can show *when* troughs
  land but cannot quantify mood (say so).

## 5. What carries over unchanged

* `population.py` / `sensitivity.py` / `personalization.py` — parameter
  names differ, machinery identical. The study layer (ingest CSV, AIMS
  analog, protocol generator) generalizes with new outcome columns.
* The **trans-masc analog of the galactagogue question** — e.g., low-dose
  T titration for nonbinary masc folks seeking partial masculinization —
  is the same "knee vs exposure" analysis shape, and there is even less
  literature. Same honesty rules apply.

## 6. Explicit non-goals / honesty rules for the wave

* No fabricated DHT/derivative modeling without citable anchors.
* Intersex populations: the original scope of this project — several
  masculinization/femininization questions (e.g., partial androgen
  insensitivity, hormone-sensitive tumor histories) need per-condition
  parameter sets, NOT one average body. Plan: condition-specific config
  packs, built only where literature anchors exist.
* Same governance everywhere: population model ≠ dosing advice; experimental
  use goes through the study layer with a clinician.

## 7. Suggested build order (next round)

1. `androgens.py` registry + T PK (esters first — the peak/trough story is
   the highest-value output)
2. T target ranges + HPG suppression anchors in `params.py`
3. `masculinization.py` (anchor on menstrual cessation; others qualitative)
4. Risk-exposure: erythrocytosis module
5. `sweep`/`compare`/`personalize` generalization pass (they're
   currently estradiol-flavored in naming)
