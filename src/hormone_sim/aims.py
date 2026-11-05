"""AIMS tracking: structured tardive-dyskinesia monitoring for a galactagogue course.

The AIMS (Abnormal Involuntary Movement Scale) is the standard structured
clinical exam for drug-induced movement disorders. It scores ten items, each
0–4 (0 none, 1 minimal, 2 mild, 3 moderate, 4 severe):

    Facial: muscles of facial expression; lips/perioral; jaw; tongue
    Extremities: upper L, upper R, lower L, lower R (chorea/athetosis)
    Trunk: neck/shoulders/hips (rocking, twisting, pelvic gyrations)
    Global: overall judgment of incapacitation

**What this module is:** a longitudinal *tracker*. It records visits, computes
totals, and raises conservative flags when a visit looks worse than the
participant's own baseline. **What it is not:** a diagnostic. TD diagnosis
belongs to the clinician; the flags exist so that worsening is impossible to
overlook, not to label anyone.

Flag logic (deliberately conservative — the failure mode is irreversible):
  * any single item ≥ 3 at any visit;
  * any item ≥ 2 points above the participant's baseline;
  * total ≥ 3 points above baseline.

Schedule encoded in :data:`MONITORING_SCHEDULE` — baseline before exposure,
every 2–4 weeks during, and at 1 and 3 months AFTER stopping (TD can emerge
or worsen post-exposure; this is the asymmetry vs lab-monitorable risks).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

AIMS_ITEMS: Dict[str, str] = {
    "facial_expression": "Muscles of facial expression (forehead, brow, smile)",
    "lips_perioral": "Lips and perioral area (puckering, smacking, pursing)",
    "jaw": "Jaw (biting, clenching, chewing, lateral movement)",
    "tongue": "Tongue (darting, rolling, protrusion)",
    "upper_left": "Upper extremities, left (chorea/athetosis)",
    "upper_right": "Upper extremities, right",
    "lower_left": "Lower extremities, left",
    "lower_right": "Lower extremities, right",
    "trunk": "Trunk — neck, shoulders, hips (rocking, twisting, gyrations)",
    "global_judgment": "Global judgment of incapacitation",
}

MONITORING_SCHEDULE = [
    ("baseline", "before first dose"),
    ("during", "every 2-4 weeks on the drug"),
    ("post-stop 1 mo", "one month after stopping"),
    ("post-stop 3 mo", "three months after stopping (TD can emerge late)"),
]


@dataclass
class AIMSVisit:
    """One AIMS administration.

    ``week`` is weeks since galactagogue start (may extend past the stop).
    Every AIMS item must be scored 0-4; use 0 for "not present" — an
    incomplete exam is recorded as such rather than silently zero-filled.
    """

    week: float
    scores: Dict[str, int]
    after_stop: bool = False
    notes: str = ""

    def __post_init__(self) -> None:
        unknown = set(self.scores) - set(AIMS_ITEMS)
        if unknown:
            raise ValueError(f"Unknown AIMS items: {sorted(unknown)}")
        missing = set(AIMS_ITEMS) - set(self.scores)
        if missing:
            raise ValueError(
                f"Incomplete AIMS exam at week {self.week}: missing "
                f"{sorted(missing)}. Score all ten items (0 if absent)."
            )
        for k, v in self.scores.items():
            if not isinstance(v, int) or not 0 <= v <= 4:
                raise ValueError(f"AIMS item {k!r} must be an int 0-4, got {v!r}")

    @classmethod
    def uniform(cls, week: float, value: int = 0, **kwargs) -> "AIMSVisit":
        """Convenience constructor: same score across all ten items."""
        return cls(week=week, scores={k: int(value) for k in AIMS_ITEMS}, **kwargs)

    def total(self) -> int:
        return int(sum(self.scores.values()))

    def worst_item(self) -> Optional[str]:
        worst = max(self.scores, key=lambda k: self.scores[k])
        return worst if self.scores[worst] > 0 else None


@dataclass
class AIMSRecord:
    """Longitudinal AIMS record for one participant."""

    baseline: AIMSVisit
    visits: List[AIMSVisit] = field(default_factory=list)
    notes: str = ""

    def __post_init__(self) -> None:
        if self.baseline.total() > 0 and not self.notes:
            # A non-zero baseline is legitimate (pre-existing movements) but
            # must be visible — flags are computed against it, not against 0.
            object.__setattr__(
                self, "notes",
                (self.notes + " | " if self.notes else "")
                + "Non-zero baseline: flags are relative to THIS baseline, "
                "and a pre-existing movement pattern warrants clinician review "
                "before starting a dopamine antagonist."
            )

    def add_visit(self, visit: AIMSVisit) -> List[str]:
        """Record a visit; return any safety flags it raises vs baseline."""
        self.visits.append(visit)
        return self.flags_for(visit)

    def flags_for(self, visit: AIMSVisit) -> List[str]:
        flags: List[str] = []
        for item, score in visit.scores.items():
            if score >= 3:
                flags.append(
                    f"week {visit.week:g}: {item} scored {score} (>=3) — "
                    "pause and reassess with the prescriber"
                )
            delta = score - self.baseline.scores.get(item, 0)
            if delta >= 2:
                flags.append(
                    f"week {visit.week:g}: {item} rose {delta} points above "
                    "baseline — schedule reassessment"
                )
        delta_total = visit.total() - self.baseline.total()
        if delta_total >= 3:
            flags.append(
                f"week {visit.week:g}: total rose {delta_total} points above "
                "baseline ("
                f"{self.baseline.total()} -> {visit.total()})"
            )
        return flags

    def report(self) -> str:
        lines = [
            "AIMS tracking record",
            f"  Items scored 0-4 x {len(AIMS_ITEMS)}; total range 0-40.",
            f"  Baseline: week {self.baseline.week:g}, "
            f"total {self.baseline.total()}",
        ]
        for v in self.visits:
            tag = " (post-stop)" if v.after_stop else ""
            lines.append(
                f"  Visit: week {v.week:g}{tag}, total {v.total()}, "
                f"worst item {v.worst_item() or 'none'}"
            )
            for f in self.flags_for(v):
                lines.append(f"    FLAG: {f}")
        lines.append("  Schedule: " + "; ".join(
            f"{name} — {when}" for name, when in MONITORING_SCHEDULE))
        lines.append(
            "  Flags are conservative safety signals, not diagnoses; "
            "interpretation belongs to the prescriber."
        )
        if self.notes:
            lines.append(f"  Notes: {self.notes}")
        return "\n".join(lines)


def parse_visit(spec: str) -> AIMSVisit:
    """Parse a CLI visit spec: 'week:all=N' or 'week:item=N,item=N[,after_stop]'.

    Examples:
        '0:all=0'
        '4:all=0'
        '8:lips_perioral=2,after_stop'
    """
    import re
    m = re.match(r"^\s*([0-9.]+)\s*:(.*)$", spec.strip())
    if not m:
        raise ValueError(f"Bad AIMS visit {spec!r}; want 'week:item=N,...'")
    week = float(m.group(1))
    after_stop = False
    scores: Dict[str, int] = {}
    uniform: Optional[int] = None
    for tok in m.group(2).split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok.lower() == "after_stop":
            after_stop = True
            continue
        k, _, v = tok.partition("=")
        k = k.strip().lower()
        if k == "all":
            uniform = int(v)
            continue
        key = k if k in AIMS_ITEMS else k.replace(" ", "_").replace("-", "_")
        if key not in AIMS_ITEMS:
            raise ValueError(f"Unknown AIMS item {k!r} in {spec!r}")
        scores[key] = int(v)
    if uniform is not None and scores:
        raise ValueError("Mix 'all=' with per-item scores in one visit, not both")
    if uniform is not None:
        return AIMSVisit.uniform(week, uniform, after_stop=after_stop)
    # CLI per-item specs zero-fill the unmentioned items: the user recording
    # 'week:lips_perioral=2' means "everything else was 0". (The AIMSVisit
    # constructor itself stays strict — programmatic callers should never
    # silently zero-fill a half-filled exam dict.)
    filled = {k: 0 for k in AIMS_ITEMS}
    filled.update(scores)
    return AIMSVisit(week=week, scores=filled, after_stop=after_stop)
