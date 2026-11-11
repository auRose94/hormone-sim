# DHT / 5α-reductase layer — references pass (AGENTS.md §5b.1)

**Date:** 2026-11-05
**Status:** references pass COMPLETE. Verdict: a **thin** DHT layer is
anchorable; a **full outcome re-driver is not**. Implementation has NOT
started — this memo is the scope guard, per the package's cite-or-skip rule
(no anchor, no model; and equally: *anchor, but only as much model as the
anchor holds*).

**Why this pass exists.** §5b.1 asked whether DHT can be modeled for the
trans-masculine drivers (voice/hair/body) and as the honest basis for 5α-RD
and PAIS intersex packs. The references pass found the anchor structure is
*narrower and differently shaped* than folk endocrinology suggests — and
found one unexpected, well-anchored trans-masculine finding that upgrades a
question raised in project discussions into a measurable proposal.

---

## 1. What 5α-reductase is (the floor)

Two isozymes convert testosterone → dihydrotestosterone (DHT): SRD5A2
(genital skin, prostate, dominant prenatally) and SRD5A1 (skin/scalp
widely). DHT is a higher-affinity androgen-receptor ligand and is **not
aromatizable** — it is the non-estrogenizable branch of the androgen
program. Precise potency multiples are assay-dependent; the package would
carry any such number as calibrated, not cited.

## 2. The natural experiment that bounds the layer (the key finding)

5α-reductase-2 deficiency (46,XY, SRD5A2 loss-of-function) is a *clean
in-vivo experiment* on what DHT uniquely does post-puberty:

* Affected individuals are undervirilized at birth (DHT-dependent external
  genital differentiation), then **virilize substantially at puberty on
  testosterone alone: muscular habitus, voice deepening, beard growth**
  (Imperato-McGinley et al., Science 1974; Thigpen et al., NEJM 1992 —
  SRD5A2 molecular basis; Imperato-McGinley & Zhu, 2002 review; StatPearls
  on 5-alpha reductase deficiency).

**Modeling consequence (load-bearing if implemented):** the masculinization
clocks in `masculinization.py` — voice, terminal hair, body composition —
should **stay testosterone-driven**. A DHT layer must NOT re-parameterize
them; the deficiency phenotype shows T carries that pubertal program with
partial compensation via the remaining isozyme. This is the same discipline
as the intersex packs: normalize against the body's own documented
physiology.

## 3. Where DHT IS the driver (the anchorable thin layer)

| Domain | Anchor | Tier |
|---|---|---|
| **Route-dependent DHT:T ratio** — transdermal gel elevates DHT via cutaneous 5α-RD (local conversion before circulation) | Wang et al. 2000 (JCEM, PMID 10720024): gel → serum DHT ~33% higher than expected; Madsen et al. 2022 (Endocrine Connections, PMC9346330) gel vs injections in trans men; Pastuszak et al. 2015 (PMC4599554) | **Cited** — supports a calibrated `dht_t_ratio` modifier per `input_kind` (gel > depot) in `androgens.py` |
| **Scalp hair (androgenic alopecia)** — DHT-mediated; inhibitor response quantifies the pathway | Finasteride ~65–70% serum DHT suppression (type-2 selective); dutasteride >90% (dual type-I/II): PMC1472914; Expert Opin Pharmacother 2023 review | **Cited** — supports a hair-note + adjunct registry candidates |
| **Prostate-analog tissue in trans men on T** — prostatic metaplasia is *documented, not hypothetical* | Xu et al. 2022 (PMID 35034167): prostatic metaplasia of the vagina in **21/21 (100%)** transmasculine vaginal specimens on testosterone; Ghent cohort: PSA expression in ~69% of metaplastic cases; rat Skene's gland acquires male-prostate phenotype on T (cited via search) | **Cited (tissue level)** — see §5 for what remains unmeasured |
| **Bottom growth / external genitalia** — genital skin is 5α-RD2-rich; DHT is the locally dominant androgen; DHT cream used for clitoral enlargement in practice | clinical-practice documents (e.g., trans-health guidelines); no dose-response anchor | **Practice-flagged** — qualitative note, no curve |
| **Sebaceous / acne** — common on T therapy | Wierckx et al. 2014 (PMID 24828032): acne, increased libido, clitoral pain on T | **Cited (incidence only)** — keep qualitative |

## 4. What the evidence says NOT to model (equally load-bearing)

1. **No serum-DHT → dermatological-outcome curve.** Wierckx et al. 2014
   (J Sex Med, "Short- and Long-Term Clinical Skin Effects of Testosterone
   Treatment in Trans Men") measured DHT and found dermatological outcomes
   (facial/body hair, Ferriman–Gallwey 0.5 → 12 median) were **not
   demonstrably related to individual serum T or DHT levels**. Serum DHT is
   a poor predictor of local tissue response — inter-individual variance
   dominates. Same honesty rule as "genetics dominate ceilings."
2. **No voice re-driver.** Deficiency phenotype (§2) virilizes voice
   without DHT.
3. **No acute "secretagogue" framing anywhere.** Prostatic-analog secretion
   is chronic androgen-trophic tissue response, not a fill-and-release tank;
   no drug creates it acutely, and emission-phase agents (α-agonists, PDE5)
   act only where the plumbing exists. The package should say so in plain
   language if this question surfaces in community materials.

## 5. The unexpected upgrade: PSA as a candidate trans-masculine biomarker

Xu 2022 establishes the *tissue* response (metaplasia in 100% of specimens).
What remains unmeasured, to the limit of this pass: **serum PSA dynamics in
trans men on testosterone**, and any relationship to functional/subjective
outcomes. This makes PSA a legitimate study-layer proposal:

* PSA is a blood biomarker (the Hct mirror, not the TD mirror: monitorable,
  plausibly dose/route-dependent, reversible on discontinuation);
* a case series could add PSA to the existing CBC draws at near-zero cost;
* the package's role: propose the measurement, model nothing until data
  exists. A candidate line for the `protocol-t` §6 measurement schedule,
  flagged as exploratory.

This resolves the trans-masculine "fluid build-and-release" question at the
level the package is honest about: the tissue is real and androgen-driven
(cited), the subjective/functional experience is community-reported
(anecdotal), and the serum biomarker is unmeasured (gap → instrument).

## 6. Intersex pack implications (§5b.6 crossover)

* **5α-RD2-deficiency pack: now anchorable.** Trajectory shape is
  documented: undervirilized at birth → pubertal virilization (voice, beard,
  muscle) on T ± residual type-1 DHT; gender-identity outcomes heterogeneous
  (Nascimento et al. 2018). Baseline overlay + anchor rating "cited
  phenotype trajectory, population-typical only."
* **PAIS:** still framework-only (receptor-sensitivity individuality stands;
  nothing in this pass changes that).
* Bonus anchor found: Wierckx 2014's Ferriman–Gallwey trajectory (median
  0.5 → 12) is a candidate real-data anchor for the `hair` outcome clock in
  `masculinization.py` — needs full-text for the time axis before use.

## 7. Implementation verdict

1. **Anchorable now (thin layer):** per-route DHT:T modifier in the
   androgen registry (calibrated), alopecia/finasteride/dutasteride note
   layer, 5α-RD2-deficiency intersex pack, protocol-t optional PSA line
   (exploratory flag).
2. **Deferred/never:** DHT-driven voice/hair re-parameterization; serum-DHT
   response curves; acute-secretion framing.
3. **Sequencing note:** SHBG (§5b.2) and estrone/first-pass (§5b.3) touch
   the same free-fraction/PK plumbing — finish all three references passes
   first, then implement as one coherent wave rather than three patches.

## Sources

* Imperato-McGinley J, et al. Steroid 5α-reductase deficiency in man. Science 1974.
* Thigpen AE, et al. Molecular basis of steroid 5α-reductase deficiency in a large Dominican kindred. NEJM 1992.
* Imperato-McGinley J, Zhu YS. Androgens and male physiology — the syndrome of 5α-reductase-2 deficiency. Mol Cell Endocrinol 2002 (PMID 12573814).
* Wang C, et al. Testosterone gel (AndroGel) pharmacokinetics. JCEM 2000 (PMID 10720024).
* Madsen MV, et al. Testosterone in men with hypogonadism and transgender males. Endocrine Connect 2022 (PMC9346330).
* Pastuszak AW, et al. Testosterone preparations for application in transgender men. 2015 (PMC4599554).
* Wierckx K, et al. Short- and long-term clinical skin effects of testosterone treatment in trans men. J Sex Med 2014 (doi:10.1111/jsm.12366); Cross-sex hormone therapy safety (PMID 24828032).
* Xu X, et al. Prostatic metaplasia of the vagina in transmasculine individuals on testosterone. 2022 (PMID 35034167).
* Zaviacic M, Ablin RJ. The female prostate and prostate-specific antigen. 2000 (PMID 10668204).
* 5α-reductase inhibitors: PMC1472914 (finasteride vs dutasteride trials); Expert Opin Pharmacother 2023 (doi:10.1080/14656566.2023.2280630).
* Nascimento R, et al. Gender identity in 5α-reductase deficiency. 2018.
* Moravek MB, et al. Effects of exogenous testosterone in trans men (review, PMC7046016).
