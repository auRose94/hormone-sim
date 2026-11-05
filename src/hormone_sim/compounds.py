"""Estrogenic compounds: a registry of real phytoestrogens plus a constructor
for *hypothetical* molecules (the "care yet to be discovered" case).

The keystone idea is **relative estrogenic activity (REA)**: every compound is
characterized by how strongly it activates the estrogen receptor relative to
17β-estradiol (= 1.0), together with its own pharmacokinetics. Two compounds
with the same plasma concentration are NOT equally estrogenic — genistein at
1 pg/mL activates far less receptor than estradiol at 1 pg/mL. REA captures
that. Combined with each compound's PK, we can put estradiol, a soy isoflavone,
a topical, and a hypothetical future drug on one common scale: **estradiol
equivalents (EE)**. See :mod:`hormone_sim.estrogenic_activity`.

REMs are *assay-dependent* and ERα/ERβ-selectivity matters; the values here are
representative literature central estimates (see docs/references.md) and are
flagged as calibrated-with-uncertainty. Treat the numbers as order-of-magnitude,
not precision.

This is pharmacology modeling for an adult researching her own care — not
dosing advice.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional


# --- reference molecular weights (g/mol) -------------------------------------
MW_ESTRADIOL = 272.38
FF_ESTRADIOL = 0.02   # ~2% free (unbound) in plasma; rest bound to SHBG/albumin


@dataclass(frozen=True)
class Compound:
    """An estrogenic compound's receptor + physical properties.

    * ``rea``            — relative estrogenic activity vs estradiol (=1.0),
                           ER-weighted intrinsic activity. Assay-dependent.
    * ``free_fraction``  — unbound (plasma) fraction available to reach tissue.
    * ``mw``             — molecular weight (g/mol), for molar↔mass conversion.
    * ``potency``        — dimensionless relative potency per unit *mass*
                           concentration, = rea * (MW_e2/MW) * (ff/ff_e2). This is
                           the number that folds a compound's plasma curve into
                           estradiol-equivalents. Cached from the other fields.
    """

    name: str
    rea: float
    free_fraction: float
    mw: float
    citation: str = ""
    calibrated: bool = True
    note: str = ""

    @property
    def potency(self) -> float:
        """Relative estrogenic potency per pg/mL vs estradiol per pg/mL."""
        return self.rea * (MW_ESTRADIOL / self.mw) * (self.free_fraction / FF_ESTRADIOL)


# --- phytoestrogen registry (real, off-the-shelf compounds) ------------------
# REA values are representative ERβ-weighted literature central estimates.
# Sources: Bovee et al. (J Agric Food Chem); Muanprasat & Chatsudthipong;
# reviews of dietary phytoestrogen endocrinology (see docs/references.md).
# free_fraction: phytoestrogens bind SHBG more weakly than E2 -> higher free frac.

PHYTOESTROGENS: Dict[str, Compound] = {
    # name:            REA      free_frac  MW      (one-line provenance)
    "genistein":  Compound("genistein",  1e-3, 0.10, 270.24,
        "Soy/red clover isoflavone; REA ~1e-3 (ERβ), Bovee et al. 2004; "
        "assay-dependent 1e-4..3e-3."),
    "daidzein":   Compound("daidzein",   1e-4, 0.10, 254.24,
        "Soy isoflavone; REA ~1e-4. Precursor to equol."),
    "equol":      Compound("equol",      5e-4, 0.10, 254.24,
        "Daidzein metabolite; ~half of populations produce it (gut flora). "
        "REA ~5e-4."),
    "glycitein":  Compound("glycitein",  1e-4, 0.10, 266.28,
        "Soy isoflavone; REA ~1e-4."),
    "formononetin": Compound("formononetin", 1e-4, 0.10, 268.26,
        "Red clover isoflavone (methoxy precursor to genistein)."),
    "biochanin_a": Compound("biochanin_a", 1e-3, 0.10, 284.26,
        "Red clover isoflavone; REA ~1e-3."),
    "coumestrol": Compound("coumestrol", 1e-2, 0.10, 286.23,
        "Coumestan; one of the most potent natural phytoestrogens, REA ~1e-2."),
    "resveratrol": Compound("resveratrol", 1e-3, 0.10, 228.24,
        "Stilbene; weak ER activity ~1e-3; also modest aromatase/AI-like effects."),
    "enterolactone": Compound("enterolactone", 3e-4, 0.10, 300.36,
        "Lignan gut-metabolite; REA ~3e-4."),
}


def get_compound(name: str) -> Compound:
    """Look up a phytoestrogen by name (raises KeyError if unknown)."""
    if name not in PHYTOESTROGENS:
        raise KeyError(
            f"Unknown compound {name!r}. Known: {sorted(PHYTOESTROGENS)}"
        )
    return PHYTOESTROGENS[name]


def hypothetical(
    name: str,
    rea: float,
    mw: float = MW_ESTRADIOL,
    free_fraction: float = FF_ESTRADIOL,
    note: str = "",
) -> Compound:
    """Construct a *hypothetical* estrogenic molecule for target-product-profile
    / inverse-design exploration ("care yet to be discovered").

    You specify the receptor activity and physical properties you're probing;
    the model then shows what PK/dose would be needed to hit a target estradiol-
    equivalent exposure. This is rational drug-design scaffolding, not a
    discovered drug.
    """
    return Compound(
        name=name, rea=rea, free_fraction=free_fraction, mw=mw,
        citation="hypothetical (user-specified target product profile)",
        calibrated=True,
        note=note or "Hypothetical compound for property-space exploration.",
    )
