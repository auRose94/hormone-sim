"""Intersex condition packs: parameter overlays for GAHT modeling.

**Why this module exists (read this first).** Intersex people — including
intersex trans people, who exist and are part of both communities — have been
systematically denied two things: bodily autonomy (non-consensual infant
"normalizing" procedures, condemned by human-rights bodies and fought by
intersex activists for decades) and **information**. Medical records
withheld, endocrinology explained to parents but not to the person it
describes. A modeling package that makes hormone physiology legible is,
for this population, not an abstraction: it is access to knowledge about
one's own body that was previously gated. This module is written in that
spirit — plain language, cited anchors, and explicit statements of what a
population model *cannot* know about a specific person.

**Language:** "intersex" (community-preferred) is used throughout; the
clinical literature also uses "differences of sex development (DSD)," which
some community members reject. Condition names below follow common usage.

**How packs work.** A pack is a *parameter overlay* on the existing models —
it changes the endogenous baselines (hormone levels without exogenous GAHT)
and endogenous production scale, and carries citations, an anchor-quality
rating, and limits. It does NOT create new biology: whatever isn't citable
isn't modeled. Every pack states what it cannot represent.

**The honesty rule, sharpened for this context:** an intersex person asking
"what will HRT do for *my* body?" has often already been given population
answers that didn't fit. These packs exist to make the model's baselines fit
the person's *documented endocrinology* — they do not pretend to replace the
clinician who knows their history.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from .androgens import BASELINE_T_NGDL, BASELINE_E2_PGML
from .params import BASELINE_E2_PGML as TF_BASELINE_E2_PGML


# Trans-feminine model baseline lives in params.py (15 pg/mL, AMAB-typical);
# androgen model baseline lives in androgens.py. Packs override per condition.

@dataclass(frozen=True)
class ConditionPack:
    """An intersex condition's parameter overlay for the GAHT models."""

    key: str
    name: str
    context: str                       # karyotype/gonadal context, plain language
    endocrinology: str                 # the cited endocrine picture
    # Parameter overlays (None = use population default):
    baseline_t_ngdl: Optional[float]   # endogenous T without exogenous GAHT
    baseline_e2_pgml: Optional[float]  # endogenous E2 without exogenous GAHT
    ovarian_production_scale: float    # 0 = no gonadal output (e.g. post-gonadectomy)
    # Trans-feminine estradiol-model baseline (may differ from the masc one):
    tf_baseline_e2_pgml: Optional[float] = None
    anchor_quality: str = "moderate"   # well-characterized | moderate | framework-only
    citations: str = ""
    limits: str = ""                   # what the model cannot represent

    def describe(self) -> str:
        lines = [
            f"{self.name} ({self.key})",
            f"  Context: {self.context}",
            f"  Endocrinology: {self.endocrinology}",
            f"  Baseline T (endogenous): "
            f"{self.baseline_t_ngdl if self.baseline_t_ngdl is not None else 'population default'} ng/dL",
            f"  Baseline E2 (endogenous): "
            f"{self.baseline_e2_pgml if self.baseline_e2_pgml is not None else 'population default'} pg/mL",
            f"  Gonadal production scale: {self.ovarian_production_scale:g}",
            f"  Anchor quality: {self.anchor_quality}",
            f"  Citations: {self.citations}",
        ]
        if self.limits:
            lines.append(f"  LIMITS: {self.limits}")
        return "\n".join(lines)


PACKS: Dict[str, ConditionPack] = {
    "cais": ConditionPack(
        key="cais",
        name="Complete Androgen Insensitivity Syndrome (CAIS)",
        context="46,XY with complete androgen-receptor insensitivity; "
        "testes may be intact or removed (gonadectomy is common, often in "
        "adolescence/early adulthood).",
        endocrinology="Testosterone in the male range but receptor-inactive; "
        "E2 arises from testicular aromatization (spontaneous breast "
        "development is typical). After gonadectomy: endogenous sex steroids "
        "drop to adrenal levels — exogenous estradiol is the long-term therapy, "
        "and estrogen-receptor function is INTACT (E2 response is normal).",
        baseline_t_ngdl=800.0,      # testicular, intact (pre-gonadectomy)
        baseline_e2_pgml=None,      # testicular aromatization while intact
        ovarian_production_scale=0.0,   # post-gonadectomy context (GAHT baseline)
        tf_baseline_e2_pgml=8.0,    # adrenal-only E2 post-gonadectomy
        anchor_quality="well-characterized",
        citations="Hughes IA et al., 'Androgen insensitivity syndrome,' "
        "Lancet 2012; Endocrine Society CPG (Grimley et al. 2017), DSD section.",
        limits="AR insensitivity is receptor-level and is NOT modeled "
        "mechanically — this pack changes endogenous baselines only. The "
        "trans-feminine estradiol model applies well post-gonadectomy (ER "
        "intact; breast development responds to E2 normally). No uterus: "
        "progesterone phases are irrelevant.",
    ),
    "pais": ConditionPack(
        key="pais",
        name="Partial Androgen Insensitivity Syndrome (PAIS)",
        context="46,XY (or 46,XX variants) with PARTIAL androgen-receptor "
        "sensitivity — the spectrum is wide and individual.",
        endocrinology="T often high (receptor resistance → compensatory "
        "gonadotropin drive); phenotype depends on the individual's actual "
        "receptor sensitivity, which varies enormously and is not predictable "
        "from population data.",
        baseline_t_ngdl=700.0,
        baseline_e2_pgml=None,
        ovarian_production_scale=0.0,
        tf_baseline_e2_pgml=15.0,
        anchor_quality="framework-only",
        citations="Hughes IA et al., Lancet 2012; consensus on DSD management "
        "(Lee PA et al., 'Global DSD update consensus,' 2016).",
        limits="FRAMEWORK ONLY: the decisive variable — individual androgen-"
        "receptor sensitivity — cannot be modeled from population parameters. "
        "Use these baselines for hormone-level reasoning; do NOT use the "
        "masculinization/femininization outcome clocks for PAIS bodies "
        "without clinician interpretation.",
    ),
    "cah_21ohd": ConditionPack(
        key="cah_21ohd",
        name="Congenital Adrenal Hyperplasia, 21-hydroxylase deficiency (CAH)",
        context="46,XX or 46,XY with an adrenal enzyme block: cortisol is "
        "low, ACTH drives adrenal androgen overproduction. Managed with "
        "glucocorticoids (and often fludrocortisone). Relevant across "
        "trans-masculine and trans-feminine care — intersex trans people "
        "exist, including on masculinizing T.",
        endocrinology="Adrenal androgen excess (androstenedione → T) raises "
        "endogenous testosterone above typical for the karyotype; "
        "glucocorticoid replacement suppresses it. Monitoring uses 17-OHP "
        "and androstenedione (not modeled here).",
        baseline_t_ngdl=150.0,      # adrenal-driven; glucocorticoid-controlled
        baseline_e2_pgml=None,
        ovarian_production_scale=1.0,   # gonads intact; adrenal adds on top
        tf_baseline_e2_pgml=None,
        anchor_quality="well-characterized",
        citations="Speiser PW et al., 'Congenital adrenal hyperplasia due to "
        "21-hydroxylase deficiency: an Endocrine Society clinical practice "
        "guideline,' 2010 (update 2018).",
        limits="The adrenal/enzyme-block structure (17-OHP, androstenedione, "
        "glucocorticoid dose-response) is NOT modeled — the pack shifts only "
        "the endogenous T baseline. Salt-wasting crises and glucocorticoid "
        "titration are entirely out of scope; endocrinology care for CAH "
        "itself continues unchanged.",
    ),
    "klinefelter": ConditionPack(
        key="klinefelter",
        name="Klinefelter syndrome (47,XXY)",
        context="47,XXY with testicular insufficiency; many are on T "
        "replacement independent of GAHT.",
        endocrinology="Primary testicular failure: low-to-low-normal T with "
        "elevated gonadotropins (LH/FSH). Trans-masculine people with "
        "Klinefelter start from a lower T baseline and titrate toward the "
        "same male-range targets as anyone else.",
        baseline_t_ngdl=250.0,      # low-male / borderline
        baseline_e2_pgml=None,
        ovarian_production_scale=1.0,
        tf_baseline_e2_pgml=None,
        anchor_quality="well-characterized",
        citations="Gravholt CH et al., 'Klinefelter syndrome: integrating "
        "genetics, neuropsychiatric, and endocrinology,' Endocr Rev 2011; "
        "Grimley et al. 2017.",
        limits="Sertoli-cell/seminal parameters, neurodevelopmental "
        "associations: out of scope. The lower baseline is the modeled "
        "difference; titration logic is unchanged.",
    ),
    "turner": ConditionPack(
        key="turner",
        name="Turner syndrome (45,X)",
        context="45,X with streak gonads (ovarian dysgenesis): no meaningful "
        "endogenous ovarian hormone production. Usually managed with estrogen "
        "replacement from adolescence; a trans-masculine person with Turner "
        "adds T to an already-replacement-based regimen.",
        endocrinology="Gonadal dysgenesis: endogenous E2 near zero without "
        "replacement; T near zero. GAHT here IS replacement — there is no "
        "endogenous production to suppress.",
        baseline_t_ngdl=5.0,
        baseline_e2_pgml=6.0,
        ovarian_production_scale=0.0,
        tf_baseline_e2_pgml=6.0,
        anchor_quality="well-characterized",
        citations="Gravholt CH et al., 'Turner syndrome: mechanisms and "
        "management,' Nat Rev Endocrinol 2019; Grimley et al. 2017.",
        limits="Congenital cardiovascular (aortic) and renal associations "
        "are major clinical considerations for ANY sex-hormone therapy in "
        "Turner — they are medical-screening matters, not modelable here; "
        "cardiology clearance is standard of care.",
    ),
}


def get_pack(key: str) -> ConditionPack:
    if key not in PACKS:
        raise KeyError(f"Unknown intersex pack {key!r}. Known: {sorted(PACKS)}")
    return PACKS[key]


def pack_summary() -> str:
    lines = [
        "Intersex condition packs (parameter overlays; see module docstring):",
        f"{'key':<12} | {'anchor':<20} | {'base T':>8} | {'base E2':>8} | name",
        "-" * 88,
    ]
    for key in sorted(PACKS):
        p = PACKS[key]
        bt = f"{p.baseline_t_ngdl:g}" if p.baseline_t_ngdl is not None else "default"
        be = (f"{p.baseline_e2_pgml:g}" if p.baseline_e2_pgml is not None
              else ("default" if p.tf_baseline_e2_pgml is None
                    else f"{p.tf_baseline_e2_pgml:g}"))
        lines.append(f"{key:<12} | {p.anchor_quality:<20} | {bt:>8} | {be:>8} | "
                     f"{p.name}")
    lines.append("")
    lines.append("Every pack states its limits in `hormone-sim intersex --name <key>`.")
    lines.append("Anchor ratings: well-characterized = cited endocrinology; "
                 "framework-only = structure without individual predictability.")
    return "\n".join(lines)
