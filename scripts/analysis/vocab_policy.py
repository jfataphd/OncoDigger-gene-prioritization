"""
vocab_policy.py -- the single, frozen, versioned collision policy of the cleaned vocabulary (Methods 3.3).

Every analysis that uses the CLEANED vocabulary takes its exclusions and suppressions from here; the analyses that deliberately use the
TRUE-RAW vocabulary (no filtering of any kind) do not call apply_suppression() and must not inherit the corrections below.

Policy history
  v3  2026-09-30  7 genes excluded entirely (GENE_EXCLUDE); 13 bare symbols suppressed at term level (the gene's full name and other
                  aliases stay matchable under the same identity).
  v4  2026-10-16  one bounded, evidence-based correction made after the exact-background rerun, on the corrected top-50 lists and the
                  flat-count control (Supplementary Methods S1, "Vocabulary correction after the exact-background rerun"):
                    bare symbols SHE, RHO, CCS, NODAL suppressed (sampled abstracts: pronoun; Rho GTPases and Spearman's rho; other
                    abbreviations such as Cronkhite-Canada syndrome and colon cancers; the adjective "nodal", 60 of 60 false);
                    aliases suppressed: RANK of TNFRSF11A (log-rank, rank-sum, Spearman rank; 14 of 14 false), nudE of NDE1 ("nude mice";
                    30 of 30 false), Lobe of AKT1S1 (anatomical lobe; 30 of 30 false).
                  Each gene keeps its approved symbol (unless suppressed above), its full approved name and its other aliases.
                  IMPACT needs no entry: the engine's own useful_gene_term() already deactivates the bare symbol.
  v5  2026-10-16  second and final step, from a whole-vocabulary screen (vocab_inventory_screen.py) of single-word aliases and ordinary-word
                  symbols, restricted to the genes that appear in any reported ranking (602 genes; 119 of 165 candidates). The decisions are
                  the rows marked "suppress" in vocab_correction_v5_table.py (rule and evidence there); they extend the v4 sets. Terms that are
                  genuine biological names (HER2, insulin, PARP, vimentin) and real gene names that are also common words with a mixed,
                  collection-dependent reading (KIT, MET, APC, POLE, CHOP, MIA, BAD, FAS) are retained and disclosed.
After v5 the vocabulary is frozen.
"""
from __future__ import annotations

POLICY_VERSION = "2026-10-16-v5"

# genes dropped entirely from every cleaned ranking (no safe alias exists); applied after the engine, so it does not change the lookups
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}

# v3/v4 decisions are listed here; v5 decisions are read from the frozen decision table (single source, no drift)
import vocab_correction_v5_table as _v5

# bare symbols suppressed as matchable single-word terms (term level)
SUPPRESS_SYMBOLS_V3 = {"GC", "HR", "HCCS", "BPIFA4P", "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON"}
SUPPRESS_SYMBOLS_V4_ADDED = {"SHE", "RHO", "CCS", "NODAL"}
SUPPRESS_SYMBOLS_V5_ADDED = {r[1] for r in _v5.ROWS if r[3] == _v5.SUP and r[0] == _v5.S} - SUPPRESS_SYMBOLS_V4_ADDED
SUPPRESS_SYMBOLS = SUPPRESS_SYMBOLS_V3 | SUPPRESS_SYMBOLS_V4_ADDED | SUPPRESS_SYMBOLS_V5_ADDED

# aliases suppressed as matchable single-word terms: gene id -> lower-case normalised terms
SUPPRESS_ALIASES: dict = {}
for _r in _v5.ROWS:
    if _r[3] == _v5.SUP and _r[0] == _v5.A:
        SUPPRESS_ALIASES.setdefault(_r[2], set()).add(_r[1])


def apply_suppression(genes_lexicon: dict, symbols=None, aliases=None, verbose: bool = True) -> dict:
    """Term-level suppression on a gene lexicon ({"single": {token: [entity, ...]}, "phrases": ...}).

    For each symbol, removes only that gene's entity from the token equal to the lower-cased symbol; for each (gene, alias) pair removes only
    that gene's entity from the alias token. Other genes that share a token keep it. Phrases are untouched. Returns a new lexicon dict."""
    symbols = SUPPRESS_SYMBOLS if symbols is None else symbols
    aliases = SUPPRESS_ALIASES if aliases is None else aliases
    single = {token: list(entities) for token, entities in genes_lexicon["single"].items()}
    changes = []

    def drop(token, gene_id, kind):
        entries = single.get(token)
        if not entries:
            return
        kept = [e for e in entries if e["id"] != gene_id]
        if len(kept) != len(entries):
            changes.append((kind, gene_id, token, len(entries), len(kept)))
            if kept:
                single[token] = kept
            else:
                del single[token]

    for symbol in sorted(symbols):
        drop(symbol.lower(), symbol, "symbol")
    for gene_id in sorted(aliases):
        for token in sorted(aliases[gene_id]):
            drop(token, gene_id, "alias")
    if verbose:
        n_alias = sum(len(v) for v in aliases.values())
        print(f"vocabulary policy {POLICY_VERSION}: {sum(c[0] == 'symbol' for c in changes)} of {len(symbols)} bare symbols and "
              f"{sum(c[0] == 'alias' for c in changes)} of {n_alias} aliases had an entry removed from the single-token lookup", flush=True)
        for kind, gene_id, token, before, after in changes:
            print(f"  {kind:6} {gene_id}: single[{token!r}] {before} -> {after} entities", flush=True)
    return {"single": single, "phrases": genes_lexicon["phrases"]}
