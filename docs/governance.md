# Governance, scope & ethics

This document defines what hormone-sim **is** and — just as importantly — what it
**is not**. It is part of the research record, not a disclaimer bolted on for
liability. Read it before publishing or adapting this work.

## 1. Scope: what this software does

hormone-sim is a **population pharmacokinetic/pharmacodynamic (PK–PD) model**. It
predicts, for a *statistical population*, how serum estradiol and testosterone
concentrations evolve over time under a described hormone regimen (route, dose,
anti-androgen). Its outputs are:

- concentration vs. time curves (population-typical or with Monte-Carlo bands),
- steady-state statistics (mean E2, fraction of a population reaching target),
- dose-response scans useful for **study design**.

It is intended for:

- **education** — understanding how route/dose/anti-androgen affect hormone levels;
- **research** — preclinical modeling, sensitivity analysis, power/sample-size
  estimation for future clinical studies;
- **study design** — exploring the dose range likely to hit guideline targets
  before a trial is run.

## 2. What this software does NOT do

- ❌ **It is not a patient simulator and does not produce individualized dosing.**
  Real gender-affirming therapy is titrated by a qualified clinician against
  *measured* serum levels, organ function, and cardiovascular/thrombotic risk.
  The population curves here will be wrong for any given individual.
- ❌ **It is not clinical advice.** Nothing here should be used to self-prescribe,
  to obtain hormones without a prescription, or to dose a real person.
- ❌ **It is not a validated medical device.** The parameters are literature-based
  and partly calibrated; the model is simplified. It has not been validated
  against a clinical cohort for the purpose of dosing decisions.

## 3. "Experimental therapies" — the governance bar

The prompt mentions *experimental* therapies. Any move from this **modeling**
toward testing in humans crosses a hard line that this software does not and
cannot clear on its own:

- **Human research requires ethics review.** Any interventional study in people
  needs approval from an institutional review board (IRB) / research ethics
  committee (REC), and a registered protocol.
- **Informed consent is mandatory** and must be genuinely informed — participants
  must understand that a regimen is experimental, the known/unknown risks, and
  that therapy would be medically supervised with monitoring.
- **Safety monitoring** (labs, cardiovascular/thrombotic risk, mental-health
  support) must be built into the protocol, not inferred from a simulation.
- **Equipoise.** An experimental regimen is only ethical where genuine
  uncertainty exists about the standard of care, and where the standard of care
  is not withheld.

This software is a **preclinical input** to that process — a way to reason about
plausible dose ranges and expected hormone exposure *before* anything touches a
person. It is not a substitute for the review, consent, and monitoring above.

## 4. Responsible-use commitments

If you build on this work:

1. **Cite the sources** in `docs/references.md` and mark which parameters are
   calibrated vs. measured.
2. **Keep the population-vs-individual distinction** explicit in any publication
   or derivative. Do not present model output as dosing guidance.
3. **Do not** use it to generate dosing instructions for named individuals, to
   bypass clinical care, or to support non-consensual or unmonitored use.
4. **Report limitations** — the model omits organ dysfunction, drug interactions
   beyond the listed anti-androgens, compliance variability, and much of HPG
   complexity (see `references.md`).
5. **Route toward care.** Where a reader needs actual therapy, the responsible
   answer is a qualified clinician / a gender-affirming care program, not a
   simulation.

## 5. Dual-use note

Hormone therapy for trans people is established, guideline-supported medical
care in most jurisdictions. This tool lowers the barrier to *understanding* the
endocrinology and to *designing* legitimate research. The same outputs could be
misread as dosing instructions; the governance above is the mitigation — keep
the population framing, keep the clinical gate (IRB + consent + monitoring), and
route individuals to care.
