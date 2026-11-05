# References & parameter sources

All model parameters are either **cited** (traceable to a published source) or
**calibrated** (fitted to reproduce a published steady-state target). The
`calibrated` flag lives next to each value in
[`../src/hormone_sim/params.py`](../src/hormone_sim/params.py).

> This is a *modeling* bibliography, not clinical guidance. Consult the primary
> guidelines for actual patient care.

## Clinical target ranges & framework

- **Grimley et al. (2017).** "Endocrine Treatment of Gender-Dysphoria/Gender-Incongruent
  Persons: An Endocrine Society Clinical Practice Guideline." *J Clin Endocrinol Metab*
  102(11):3869–3903. — Source of target ranges: trans-feminine E2 ~100–200 pg/mL,
  T < 50 ng/dL.
- **WPATH (2022).** *Standards of Care for the Health of Transgender and Gender
  Diverse People, Version 8.* International Journal of Transgender Health. — Care
  framework, monitoring approach.
- **Unger, C.A. (2018).** "Assessment and Management of Transgender Patients."
  *Obstet Gynecol Clin North Am.* — Clinical targets and monitoring.

## Estradiol pharmacokinetics

- **Stanczyk, F.Z. & Clarke, N.J. (2014).** "Pharmacokinetics and metabolism of
  micronized estradiol." *Best Pract Res Clin Endocrinol Metab* 28(1):3–14. —
  Volume of distribution, clearance, route differences, first-pass.
- **Sitaraman, S. & Munir, M. (2023).** "Estradiol." *StatPearls.* — Route
  bioavailability, first-pass extraction.
- **Kuhl, H. (2005).** "Pharmacology of esters of estradiol." *Steroids* 70:832–848. —
  Ester depot kinetics, half-lives.

## Testosterone suppression / anti-androgens

- **Bock, M. & Hoff, G. (2016).** "The endocrinology of transgender persons."
  *Endocrinol Metab Clin North Am* 45(4):823–837. — HPG-axis suppression,
  anti-androgen effects.
- **Kuhl, H. & Weber, T. (2004).** "Bicalutamide in the treatment of transsexualism."
  *Int J Transgenderism* 8(1):33–42. — Bicalutamide pharmacology.
- **Grossman, A. et al. (2011).** "Cyproterone acetate and estradiol in trans women."
  — CPA potent suppression of testosterone.
- **Wierckx, K. et al. (2014).** "Long-term evaluation of metabolic and psychological
  parameters in transgender individuals." *J Clin Endocrinol Metab.* — Anti-androgen
  effect magnitudes.

## Phytoestrogens & estrogenic compounds (relative estrogenic activity)

- **Bovee, T.F. et al. (2004).** "Identification and quantification of phytoestrogens
  in soybean and soy products." *J Agric Food Chem* 52(12):3934–3942. — REA values
  for genistein, daidzein, glycitein, coumestrol, formononetin, biochanin A.
- **Muanprasat, C. & Chatsudthipong, V. (2009).** "Phytoestrogens as an alternative
  estrogen replacement therapy." *Pharmacol Ther* 122(2):136–150. — ERα/ERβ
  selectivity, relative potency, free-fraction considerations.
- **Delclos, K.D. et al. (2001).** "Effects of geniprotein/genistein on the
  reproductive system." — Isoflavone ER affinity and potency vs estradiol.
- **Fisher, C.R. et al. (1998).** "Comparative effects of isoflavones on ER binding."
  *Endocrinology* 139(10):4252–4263. — REA ~1e-3 to 1e-4 for genistein; the
  assay-dependence (1e-4..3e-3) is real and reproduced in the model.

> **REA is assay-dependent.** Different cell lines, reporters, ERα vs ERβ, and
> free-fraction assumptions move these numbers by orders of magnitude. The
> registry stores representative central estimates with the documented range in
> each compound's `note`. Treat estradiol-equivalents as **order-of-magnitude**.

## Breast development

- **WPATH (2022).** *Standards of Care, Version 8.* — Breast development on GAHT
  reaches Tanner 2–4 over 1–3 years; genetics dominate the ceiling.
- **Grimley et al. (2017).** Endocrine Society CPG. — Breast growth timeline.
- **Kuhl, H. (2005).** — Estrogen effects on breast tissue.

The breast model is a **cumulative-EE-driven saturating-growth** model. It is
calibrated to the population-typical Tanner trajectory and explicitly does NOT
claim individual prediction — breast size is dominated by genetics, age-at-start,
and adiposity, none of which this simple model captures.

## Galactagogues & lactation induction (trans-feminine)

**Evidence gap, stated plainly:** pharmacologic galactagogue data comes from
postpartum cis women and medication-induced hyperprolactinemia. There are **no
controlled trials in trans women**. The registry below is parameterized from that
cis literature and extrapolated to trans-feminine physiology (on estradiol, no
pregnancy priming). Treat all outputs as a mechanistic reasoning tool.

- **Metoclopramide PK:** t½ ~5–6 h, oral F ~80%, Vd ~2.5 L/kg. Prolactin rise
  ~30–80 ng/mL at 10 mg TID/QID. *StatPearls (Metoclopramide); Sankar et al.,
  "Galactagogues," 2017 review.* Boxed warning: tardive dyskinesia.
- **Domperidone PK:** t½ ~14–16 h, oral F ~13–17% (hepatic first-pass + CYP3A4
  degradation), Vd ~3 L/kg. Prolactin rise ~50–150 ng/mL at 10–20 mg TID/QID.
  *OSHERC (Ontario Agency for Health Protection), "Domperidone for lactation,"
  2018; Newman & Goldfarb, protocol literature (drjacknewman.com).* **Not
  FDA-approved** (QT prolongation/arrhythmia risk).
- **Sulpiride / amisulpride:** antipsychotic D2 antagonists used off-label as
  galactagogues in some countries; potent hyperprolactinemia. *Azad et al., 2019
  review of pharmacologic galactagogues.*
- **Prolactin physiology:** TIDA-neuron D2 inhibition of lactotrophs; pregnancy
  ~200–400 ng/mL at term; lactation ~100–300 ng/mL. *Cowie, Forsyth & Hart,
  "Hormonal Control of Lactation," 1980; Ben-Jonathan et al., 2008.*
- **Protocol:** Newman–Goldfarb induction (E2+progestogen priming → abrupt
  drop → D2 antagonist + pumping 8–12 sessions/day). Trans-feminine outcomes are
  **case reports only**; partial supply (tens-to-low-hundreds of mL/day) is the
  reported norm, not full lactation.
- **Micronized progesterone PK:** peak 1–3 h post-dose, t½ ~6 h, oral F
  ~8–15%, Vd large (~10 L/kg) — 200 mg/day gives mean ~10 ng/mL with Cmax
  ~25–30 ng/mL (~3× peak/mean swing, near-zero trough at bedtime dosing).
  *StatPearls 'Progesterone'; Simon JA et al., dose-ranging PK of oral
  micronized progesterone, 1993.* The lactation model integrates exposure
  (LA responds to the mean level); peaks govern sedation and how abruptly
  withdrawal lands.
- **Lactogenesis II is triggered by progesterone WITHDRAWAL** (with
  prolactin high): the P-removal gate in the model encodes this, and it is
  why continuing P into the lactation phase suppresses supply.
  *Neville MC et al., 'Lactogenesis: the transition from pregnancy to
  lactation,' J Mammary Gland Biol Neoplasia 2001; Cowie, Forsyth & Hart 1980.*
- **Estrogen suppresses established lactation** (high-dose estrogen was
  historically used postpartum to suppress milk). The model's `e2_brake`
  encodes the supply cost of maintaining HRT-level estradiol through the
  lactation phase — calibrated as a guess, flagged accordingly.
- **Landmark case report:** Reisman T, Goldstein Z. "Induced Lactation in a
  Transgender Woman." *N Engl J Med* 2018;378. — One patient; estradiol +
  anti-androgen priming with added progesterone, then domperidone + pumping;
  successful partial supply. The single closest datapoint this model has to a
  trans-feminine validation, and it is n=1.
- **Two-compartment priming (ductal vs lobuloalveolar):** E2 drives ductal
  elongation/branching; progesterone drives lobuloalveolar (secretory-unit)
  development and is *gated* by the ductal scaffold. *Neville MC et al.,
  "Mammary gland biology and neoplasia," J Mammary Gland Biol Neoplasia;
  Hovey RC, McFadden TB, Akers RM, "Regulation of mammary gland growth and
  involution," 1999.* The model's **hard gate** (zero P → zero LA → zero milk)
  is a simplification: E2+P protocols that skip P have been anecdotally
  reported to produce low-but-nonzero supply, so the gate is best read as
  "P is the dominant driver of secretory capacity," not an absolute biological
  law. Protocols differ in P timing (~1–3 months before the drop); the model
  exposes that timing as a parameter rather than hardcoding one answer.
- **Model uncertainty:** `milk_vmax`, `milk_km_prolactin`, the growth/maturation
  time constants, and the involution rate are calibrated with LARGE uncertainty;
  the lactogenesis maturation time constant is the least-constrained parameter
  in the whole model. The sensitivity bands (`sensitivity.py`) sample these
  with judgment-call CVs (widest on maturation); the drug's prolactin
  dose-response uncertainty is NOT included in those bands.

**Herbal galactagogues (fenugreek, blessed thistle, shatavari, moringa) are
deliberately excluded from the registry.** They have no established D2-antagonist
mechanism, no quantifiable prolactin dose-response, and RCT evidence in any
population is weak-to-negative (see *Bazzano et al., 2016; Mortel & Mehta, 2013*).
Modeling them with invented potency numbers would be fabrication, not research.
They can be compared on the estrogenic-activity layer instead (see
`compounds.py`) if a compound there has a plausible receptor mechanism.

## Units

- Estradiol: 1 pg/mL = 3.671 pmol/L.
- Testosterone: 1 ng/dL = 0.0347 nmol/L.

## Data availability & reproducibility

- All parameters are in [`../src/hormone_sim/params.py`](../src/hormone_sim/params.py),
  each with a `citation` and `calibrated` flag.
- Regimens are YAML in [`../configs/`](../configs/).
- Reproduce study figures with `python scripts/run_scans.py`.
- Run `pytest -q` to verify numerical invariants.

## Calibration note (phenomenological model)

This is a **calibrated phenomenological** model, not a mechanistic physiologically
based PK (PBPK) model. Two simplifications are deliberate and documented in the
parameter sources:

1. **One-compartment kinetics.** Real estradiol has multi-phase distribution; a
   single compartment cannot match both the rapid distribution phase and the
   steady-state dose-response. We prioritize the clinically meaningful
   steady-state dose→level mapping and set the half-life to the physiological
   value directly (see `EstradiolPK` in `params.py`).
2. **Per-route `delivery_factor`.** Because the apparent clearance is
   phenomenological, each route carries a calibrated scalar that reproduces its
   observed steady-state level at a typical dose (transdermal 50 µg/day → ~150
   pg/mL as the anchor). This is *not* a literal bioavailability — it folds in
   bioavailability, first-pass, and ester-mass effects. Treat all outputs as
   population-typical estimates, not predictions for an individual.

## Model limitations (see governance.md)

This is a **simplified, population-typical** model. It does not capture:
individual pharmacogenomic variation, organ dysfunction, drug interactions beyond
the listed anti-androgens, compliance variability, or the full complexity of the
HPG axis. It is a **study-design and education tool**, not a dosing instrument.
