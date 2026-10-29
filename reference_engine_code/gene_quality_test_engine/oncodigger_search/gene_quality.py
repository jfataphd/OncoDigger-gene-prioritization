"""Auditable gene vocabulary policies shared by the app and offline analyses."""
from pathlib import Path

import pandas as pd

from .lexicons import add_entity, build_lookup, norm_text, useful_gene_term

POLICY_VERSION = "2026-09-27-v1"
REVIEW_MULTIPLIER = 0.65
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
EXCLUSIONS = {
    "MICE": "Collides with mice in animal-model reports",
    "DCR": "Collides with disease control rate",
    "PC": "Collides with prostate cancer / principal component",
    "C2": "Collides with structural and other abbreviations",
    "FBN1": "Existing app QC policy: reported false-positive cervical-cancer matches",
    "SGCG": "Existing app QC policy: reported non-gene alias collisions",
    "WDHD1": "Existing app QC policy: reported query-independent distributional artifact",
}
SAFE_ALIASES = {"GC": "VDBP", "HR": "KDM3D", "HCCS": "CCHL", "BPIFA4P": "LATH"}


def gene_lookups(root: Path) -> tuple[dict, dict, list[dict]]:
    """Raw uses all supplied HGNC terms; Vetted applies explicit, logged QC.

    Neither mode changes document retrieval or disables correctness fixes.
    Canonical HGNC identity is retained even when a safer alias is displayed.
    """
    frame = pd.read_csv(root / "Human Genes.csv").fillna("")
    approved = {norm_text(s): str(s).upper() for s in frame["symbol2"] if s}
    raw, vetted, audit = {}, {}, []

    def note(symbol, term, action, reason):
        audit.append({"HGNC Symbol": symbol, "Term": term, "Action": action, "Reason": reason})

    for row in frame.to_dict("records"):
        original = str(row["symbol2"] or row["symbol"]).upper()
        symbol = GENE_REMAP.get(original, original)
        name = str(row["Approved name"])
        aliases = [t.strip() for t in str(row["Alias symbols"]).split(",") if t.strip()]
        aliases.extend(alias for alias, target in GENE_REMAP.items() if target == symbol)
        terms = sorted({norm_text(t) for t in [original, symbol, name, *aliases] if norm_text(t)})
        # Merge obsolete identity rows into the approved symbol before scoring.
        if symbol in raw:
            terms = sorted(set(terms) | set(raw[symbol]["terms"]))
        add_entity(raw, symbol, symbol, "Gene", terms, source_id=symbol, description=name)
        if symbol in EXCLUSIONS:
            note(symbol, symbol, "Excluded gene", EXCLUSIONS[symbol])
            continue
        display = SAFE_ALIASES.get(symbol, symbol)
        accepted = []
        if display != symbol:
            note(symbol, display, "Safer alias", f"Display {display}; suppress ambiguous primary symbol {symbol}")
        for term in terms:
            owner = approved.get(term)
            if owner and GENE_REMAP.get(owner, owner) != symbol:
                note(symbol, term, "Removed alias", f"Collides with approved symbol {owner}")
            elif display != symbol and term == norm_text(symbol):
                note(symbol, term, "Removed term", "Ambiguous primary symbol; use safer alias or full name")
            elif (term in {norm_text(a) for a, target in GENE_REMAP.items() if target == symbol}
                  or useful_gene_term(term, display)):
                accepted.append(term)
            else:
                note(symbol, term, "Removed term", "Short or common-word alias filter")
        if symbol in vetted:
            accepted = sorted(set(accepted) | set(vetted[symbol]["terms"]))
        add_entity(vetted, symbol, display, "Gene", accepted, source_id=symbol, description=name)
    for alias, symbol in GENE_REMAP.items():
        note(symbol, alias, "Canonical alias", f"Score {alias} mentions under {symbol}")
    return build_lookup(raw), build_lookup(vetted), audit
