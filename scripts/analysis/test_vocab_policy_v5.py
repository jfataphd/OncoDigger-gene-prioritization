"""
test_vocab_policy_v5.py -- preflight test of the v5 vocabulary policy at the matcher level (Supplementary Methods S1).

Replays the engine's own per-paper matching block (oncodigger/search/discovery.py, associated_entity_tables: token lookup in lookup["single"],
phrase lookup in lookup["phrases"]) on constructed texts, with the raw gene lexicon and with the v5-suppressed lexicon, and checks that
(1) every suppressed term stops matching its gene, (2) the same gene is still matched through its approved symbol or full name,
(3) a different gene that shares the suppressed token is unaffected, (4) retained terms (including the eight retain-collection-dependent ones)
still match, (5) the raw lexicon is untouched by apply_suppression.
Run: ONCODIGGER_ROOT=<frozen corpus root> python -X utf8 scripts/analysis/test_vocab_policy_v5.py
"""
import os, sys
from pathlib import Path
ONC = Path(os.environ["ONCODIGGER_ROOT"])
sys.path.insert(0, str(ONC / "src"))
import vocab_policy as vp
import vocab_correction_v5_table as t5
from oncodigger.search.lexicons import load_lexicons, norm_text
from oncodigger.search.discovery import full_document_text, document_tokens

raw = load_lexicons(ONC)["Genes"]
raw_single_snapshot = {k: [e["id"] for e in v] for k, v in raw["single"].items()}
sup = vp.apply_suppression(raw, verbose=False)
assert {k: [e["id"] for e in v] for k, v in raw["single"].items()} == raw_single_snapshot, "apply_suppression modified the raw lexicon"

def genes(lex, text):
    doc = {"title": text, "abstract": ""}
    padded = f" {full_document_text(doc)} "
    out = set()
    for tok in set(document_tokens(doc)):
        out |= {e["id"] for e in lex["single"].get(tok, [])}
        for phrase, e in lex["phrases"].get(tok, []):
            if f" {phrase} " in padded:
                out.add(e["id"])
    return out

fails = []
def check(name, cond):
    print(("ok   " if cond else "FAIL ") + name)
    if not cond: fails.append(name)

# (1)+(3): each suppressed (gene, term) stops matching; other genes sharing the token keep it
for kind, term, gene, dec, *_ in t5.ROWS:
    if dec != t5.SUP: continue
    text = f"the {term.lower()} was reported here"
    check(f"suppressed {kind} {term!r} -> {gene}: raw matches, v5 does not", gene in genes(raw, text) and gene not in genes(sup, text))
shared = {"ca125": ("MUC16", "NBR1"), "cd25": ("IL2RA", "ISG20"), "nrf2": ("NFE2L2", "GABPA")}
for tok, (keep, drop) in shared.items():
    g = genes(sup, f"the {tok} marker")
    check(f"shared token {tok!r}: {keep} kept, {drop} dropped", keep in g and drop not in g)
# (2): the same gene is still reachable through a remaining route (another single-word name, or a full-name phrase)
def routes(lex, gene):
    singles = sorted(t for t, v in lex["single"].items() if any(e["id"] == gene for e in v))
    phrases = sorted(p for l in lex["phrases"].values() for p, e in l if e["id"] == gene)
    return singles, phrases
suppressed_genes = sorted({r[2] for r in t5.ROWS if r[3] == t5.SUP})
for gene in suppressed_genes:
    singles, phrases = routes(sup, gene)
    probe = [f"the {x} was reported" for x in singles[:2]] + [f"the {p} was reported" for p in phrases[:2]]
    ok = any(gene in genes(sup, t) for t in probe)
    check(f"{gene}: still matched after suppression via {len(singles)} single / {len(phrases)} phrase routes", ok or (not singles and not phrases and False))
# (4): retained terms still match
for kind, term, gene, dec, *_ in t5.ROWS:
    if dec in (t5.RET, t5.RCD):
        check(f"retained {term!r} -> {gene} ({dec})", gene in genes(sup, f"the {term.lower()} was reported"))
# retain-collection-dependent is documentation only: none may appear in the suppression sets
rcd_syms = {r[1] for r in t5.ROWS if r[3] == t5.RCD and r[0] == t5.S}
rcd_alias = {(r[2], r[1]) for r in t5.ROWS if r[3] == t5.RCD and r[0] == t5.A}
check("retain-collection-dependent symbols not in SUPPRESS_SYMBOLS", not (rcd_syms & vp.SUPPRESS_SYMBOLS))
check("retain-collection-dependent aliases not in SUPPRESS_ALIASES", not any(t in vp.SUPPRESS_ALIASES.get(g, ()) for g, t in rcd_alias))
print("\nFAILURES:", fails if fails else "none")
sys.exit(1 if fails else 0)
