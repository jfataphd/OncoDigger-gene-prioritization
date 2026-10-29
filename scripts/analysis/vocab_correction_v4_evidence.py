"""
vocab_correction_v4_evidence.py -- evidence record for the v4 vocabulary correction (vocab_policy.py; Supplementary Methods S1).

For every term that was corrected, and for the audited genes that were left alone, this regenerates the seeded random sample of papers in
which the term occurs (random.seed(42), sampled from the token's postings in the collection where the gene was flagged) and records the
collection, the number of papers with the token, the sample PMIDs, and the reading of that sample (n_genuine = abstracts in which the term
refers to the gene; the rest use an unrelated meaning). The readings were made by a language-model agent that had seen the manuscript, so
they are not independent validation. No abstract text is stored; the PMID list lets a reader re-read the same abstracts.

Output: data/derived/vocab_rescoring/vocab_correction_v4_evidence.csv
Run: ONCODIGGER_ROOT=<frozen corpus root> python -X utf8 scripts/analysis/vocab_correction_v4_evidence.py
"""
from __future__ import annotations

import csv
import os
import random
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORPUS = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_processed"
OUT = ROOT / "data" / "derived" / "vocab_rescoring" / "vocab_correction_v4_evidence.csv"

# term, gene, kind, collection, sample size, genuine mentions read, meaning of the non-gene uses, decision
ROWS = [
    ("SHE", "SHE", "symbol", "thyroid", 8, 0, "the pronoun (case reports: 'she underwent ...')", "suppress bare symbol"),
    ("RHO", "RHO", "symbol", "stomach", 8, 0, "Rho GTPases and Spearman's rho", "suppress bare symbol"),
    ("RHO", "RHO", "symbol", "stomach", 30, 0, "Rho GTPases (RhoA, Rac1, Cdc42, RhoGDI2, ARHGAP) and Spearman's rho", "suppress bare symbol"),
    ("CCS", "CCS", "symbol", "colorectal", 8, 0, "colorectal cancer screening, Charlson and Cleveland Clinic scores, 'CCs' (colon cancers)", "suppress bare symbol"),
    ("CCS", "CCS", "symbol", "colorectal", 30, 3, "Cronkhite-Canada syndrome, colon cancers (CCs), scoring systems, a stapler, screening; 3 of 30 were the copper chaperone gene", "suppress bare symbol"),
    ("NODAL", "NODAL", "symbol", "thyroid", 30, 0, "the adjective 'nodal' (nodal metastasis, nodal dissection)", "suppress bare symbol"),
    ("NODAL", "NODAL", "symbol", "lymphoma", 30, 0, "the adjective 'nodal' (nodal and extranodal lymphoma)", "suppress bare symbol"),
    ("rank", "TNFRSF11A", "alias", "breast", 14, 0, "log-rank, rank-sum, Spearman rank, ranking of interventions", "suppress alias"),
    ("nude", "NDE1", "alias", "lung", 30, 0, "nude mice and nude rats (xenograft models)", "suppress alias"),
    ("lobe", "AKT1S1", "alias", "liver", 30, 0, "anatomical lobe (hepatic lobe, caudate lobe, lung lobe)", "suppress alias"),
    # audited and left unchanged (flagged by the collision screen)
    ("ID3", "ID3", "symbol", "lymphoma", 10, 10, "", "no change"),
    ("SMO", "SMO", "symbol", "brain", 10, 10, "", "no change"),
    ("PTPRD", "PTPRD", "symbol", "melanoma", 7, 7, "", "no change"),
    ("AXIN2", "AXIN2", "symbol", "endometrial", 5, 5, "", "no change"),
    ("TRAF2", "TRAF2", "symbol", "myeloma", 8, 8, "", "no change"),
]


def main() -> None:
    out = []
    for term, gene, kind, coll, n, genuine, meaning, decision in ROWS:
        con = sqlite3.connect(f"file:{(CORPUS / f'{coll}_cancer_inclusive_all' / 'lexical_index.db').as_posix()}?mode=ro&immutable=1", uri=True)
        pm = [str(r[0]) for r in con.execute("select distinct pmid from postings where token = ?", (term.lower(),))]
        random.seed(42)
        sample = random.sample(pm, min(n, len(pm)))
        out.append({"term": term, "gene": gene, "kind": kind, "collection": coll, "papers_with_token": len(pm), "sample_n": len(sample),
                    "genuine_gene_mentions": genuine, "other_meaning": len(sample) - genuine, "other_meaning_is": meaning,
                    "decision": decision, "sample_pmids": ";".join(sorted(sample))})
        print(f"{term:6} {gene:10} {coll:11} papers {len(pm):5}  sample {len(sample):2}  genuine {genuine:2}  -> {decision}", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
