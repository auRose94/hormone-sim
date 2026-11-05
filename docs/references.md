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
