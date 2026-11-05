# AGENTS.md — hormone-sim

**For agents and human collaborators working in this repository.** Read this
before writing code. It records what exists, the load-bearing modeling
decisions (and why they are the way they are), the non-negotiable rules, and
where the research goes next. Pair with: [`README.md`](README.md) (overview),
[`docs/references.md`](docs/references.md) (parameter provenance),
[`docs/governance.md`](docs/governance.md) (scope & ethics),
[`docs/roadmap-trans-masc.md`](docs/roadmap-trans-masc.md) (masc wave status).

---

## 1. What this project is

A reproducible **population pharmacokinetic/pharmacodynamic (PK–PD) modeling
package** for gender-affirming hormone therapy, in both directions —
trans-feminine (estradiol, estrogenic compounds, breast development,
galactagogues/lactation) and trans-masculine (testosterone, masculinization,
erythrocytosis) — plus **intersex condition packs** that adapt the models to
documented intersex endocrinology, and a **study layer** (CSV ingest,
AIMS monitoring, personal inference, IRB-style protocol generation) for
clinician-supervised titration research.

It exists because there is **not enough research** on how GAHT works in
trans and intersex bodies, and because the people who need that research
have historically been denied both autonomy and information. The package's
framing is **comfort and access** — reducing internal grief by making one's
own endocrinology legible — **not "passing,"** which is other people's
opinion and not a modeling target.

## 2. Non-negotiable rules (agents: this section is binding)

1. **No fabricated parameters.** Every constant needs a citation or a
   `calibrated=True` flag with its reasoning recorded next to it. If no
   honest anchor exists, *exclude and say so in code* — precedent: herbal
   galactagogues (no D2 mechanism), oral testosterone undecanoate
   (food-dependent variable absorption), PAIS (framework-only rating).
   Inventing numbers is fabrication dressed as research.
2. **Population ≠ individual.** Model output is population-typical
   estimation, never dosing advice. The distinction must survive into every
   derivative, plot label, and report string.
3. **Bands over points.** Where uncertainty exists (everywhere), the band,
   sweep, or sensitivity envelope is the honest product; single lines are
   illustrations. Parameter CVs are judgment calls and must be labeled as
   such in output (`sensitivity.py` does this — follow the pattern).
4. **Language.** "Intersex" (community-preferred); outcomes are framed as
   comfort/access/expectation-setting, never as "passing." Genetics dominate
   outcome ceilings — say so in the output itself.
5. **Tests before claims.** `pytest -q` must pass (and new invariants
   belong in `tests/`) before any result is reported. The tests have caught
   real modeling bugs, not just typos — trust them.
6. **State limits in output.** Every model report carries its own caveat
   lines (evidence gaps, extrapolations, "not medical advice"). Preserve
   them; extend the pattern to new modules.
7. **Privacy.** Participant logs are identifiable health data: `examples/`
   contains only synthetic records; never commit real ones.

## 3. Module map

| Domain | Modules | One-line purpose |
|---|---|---|
| Trans-fem core | `params.py`, `pk.py`, `pd.py`, `regimens.py`, `sim.py`, `population.py`, `plot.py` | Estradiol 1-cmt PK (4 routes, calibrated per-route `delivery_factor`), HPG T-suppression, Monte-Carlo bands |
| Estrogenic activity | `compounds.py`, `estrogenic_activity.py`, `breast_development.py` | Phytoestrogen registry (REA), estradiol-equivalents, inverse design; EE-driven breast model |
| Galactagogues / lactation | `galactagogues.py`, `progestogens.py`, `lactation.py`, `comparison.py`, `risk_exposure.py` | D2-antagonist registry; Newman–Goldfarb protocol (two-compartment priming, P-withdrawal gate, taper, E2 brake); equieffective compare; low-dose sweep; TD exposure framing |
| Trans-masc | `androgens.py`, `masculinization.py`, `erythrocytosis.py` | T ester/gel PK (peak–trough quantified), reversed HPG suppression, outcome clocks (cessation anchor ~5 mo), CBC-monitorable risk framing |
| Intersex | `intersex.py` | Five cited condition packs (CAIS, PAIS*, CAH, Klinefelter, Turner) as baseline overlays with anchor ratings + limits (*PAIS framework-only) |
| Study layer | `measurement.py`, `personalization.py`, `personalization_androgen.py`, `ingest.py`, `aims.py`, `protocol.py`, `aggregation.py` | Domain-generic measurement-model engine (importance sampling, predictive tables, population sweeps); per-domain adapters (lactation + T logs); fail-closed CSV ingest; AIMS TD tracking; protocol skeletons for both domains (live expected values); de-identified cohort summaries |
| CLI | `cli.py` | `run scan pop compound hypothetical breast galactagogues lactation compare sweep risk personalize ingest aims protocol androgens masc hct intersex ingest-t protocol-t sweep-t cohort` |

Version at writing: **0.10.0**, **196 passing tests** across
`tests/test_hormone_sim.py`, `test_estrogenic.py`, `test_galactagogues.py`,
`test_lactation_v2.py`, `test_lowdose_risk.py`, `test_v3_lactation.py`,
`test_study_tools.py`, `test_masc.py`, `test_intersex_t.py`,
`test_t_figure.py`, `test_aggregation.py`, `test_measurement.py`.

## 4. Load-bearing modeling decisions (do not silently change)

These encode hard-won calibration lessons; several were *caught by tests*:

- **Per-route `delivery_factor` is calibrated, not literal bioavailability.**
  A one-compartment model cannot simultaneously match real estradiol's
  kinetics and its observed dose→level mapping; the transdermal 50 µg/day →
  ~150 pg/mL anchor defines apparent clearance, and other routes get
  calibrated scalars. Documented in `EstradiolPK`'s docstring.
- **Half-life and dose→level mapping are decoupled** (`half_life_h_val` +
  `steady_state_factor`) for the same reason.
- **HPG suppression is normalized against each body's OWN baseline**
  (`simulate_androgen`): zero-dose → *that body's* documented baseline for
  any pack. Un-normalized, CAH's adrenal baseline collapsed under the
  model's own feedback — the test suite caught it; keep it that way.
- **Priming and maturation are separate timescales.** Pumping *maintains*
  glandular priming but never *rebuilds* it (a test enforces this); lactosis
  II waits for progesterone **withdrawal** (`p_withdrawal_frac`), which is
  why "the drop" matters.
- **The maturation prolactin threshold is sampled** — otherwise low galactagogue
  doses show a fake cliff instead of the real (bimodal, heterogeneous) zone.
- **Estradiol-equivalents are order-of-magnitude**: REA values are
  assay-dependent by orders of magnitude; never present EE as precision.
- **The TD vs erythrocytosis asymmetry is deliberate**: TD has no biomarker
  (AIMS exam, cumulative, potentially irreversible) — `risk_exposure.py`
  outputs exposure *tiers*, never probabilities. Hct is the mirror case:
  biomarker exists, dose/route-dependent, reversible. Keep the two framings
  distinct.
- **Personal inference is importance sampling** with ESS diagnostics —
  a few observations leave posterior ≈ prior, and the tool says so.
- **The study layer's engine is domain-generic (`measurement.py`); domains
  are descriptors.** Personalization, predictive tables, and population
  sweeps run once, over a `MeasurementDomain` adapter (sample body /
  log-likelihood / curve names). Determinism contract: domains draw RNG
  *only* in `sample_body` — never in `log_likelihood` — so results stay
  reproducible and the behavior-preservation tests
  (`tests/test_measurement.py`, which re-implement the pre-abstraction
  loops and require bit-for-bit agreement) keep new domains honest.

## 5. Where the research continues (open items at time of writing)

### 5a. Engineering completions — DONE (v0.10.0; no new science was needed)

1. ✅ **Measurement-model abstraction**: the study layer's importance-
   sampling engine, posterior-predictive tables, and population sweeps are
   domain-generic (`measurement.py`); lactation and T domains are adapters
   (`LactationDomain`, `AndrogenDomain`); behavior preserved bit-for-bit
   against reference-loop tests. `sweep-t` is the first T-side sweep.
   (`compare` remains lactation-presenting by design — equieffective-
   prolactin matching is a trans-fem compound-selection question.)
2. ✅ **Masc protocol generator**: `protocol-t` / `generate_androgen_protocol`
   — T titration + CBC schedule, live trough/mean/peak + Hct-plateau
   expected values, gel-vs-depot language.
3. ✅ **Personalized figure for T logs**: `ingest-t --out` writes the
   posterior-vs-prior T/Hct figure (`plot_personal_androgen`).
4. ✅ **Aggregation view**: `cohort` / `aggregation.py` — de-identified
   multi-participant summaries from a study-ID manifest, feeding the
   protocols' statistical plans.

### 5b. Scientific expansions (need literature anchors first — cite or skip)

1. **DHT layer** (5α-reductase): voice/hair/body drivers and the honest
   basis for 5α-RD and PAIS packs. No anchor, no model — current stance.
2. **SHBG / free-hormone layer**: explicit free-vs-bound compartment
   (would refine REA/potency and T trough behavior).
3. **Estrone / oral-route metabolites**: oral E2's estrone contribution and
   VTE-relevant first-pass story (currently folded into `delivery_factor`).
4. **Prolactin dose-response in trans-feminine physiology**: unmeasured;
   the lactation layer's biggest single unknown — a study could measure it.
5. **Peak/trough symptom correlation**: patient-reported patterns are
   real but unquantified; would need EMA-style data collection (a study
   instrument, not a literature parameter).
6. **More intersex packs** (5α-RD, Swyer, MRKH, hypogonadotropic
   hypogonadism, Turner-on-T specifics): only as citations permit.

### 5c. Study/community layer

1. Run the layer for real: the protocol skeleton +
   briefing (`docs/briefing-metoclopramide.md`) are designed for a
   trans-inclusive clinic to adopt; AIMS/ingest are the instruments.
2. **Zenodo DOI + docs site** when publishing; cite the Endocrine Society
   CPG and WPATH SOC-8 alongside the repo.
3. Translations and community review of the plain-language caveats — the
   wording was written to be read by the people it describes.

## 6. Governance (short form — full text in `docs/governance.md`)

Research/education tool. Not a patient simulator, not dosing advice, not a
medical device. Human research requires IRB/REC approval and informed
consent — the study layer exists to *support* that path, never to bypass it.
Population framing + clinical gate + monitoring: keep all three.
