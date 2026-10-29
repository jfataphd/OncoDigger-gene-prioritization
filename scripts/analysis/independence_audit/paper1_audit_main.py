"""Paper 1 / IntOGen independence audit: 190-pair forensic evidence tables and benchmark recalculation.

READ-ONLY on Paper 1 (OncoDigger-mutation-paper1 working tree) and on the IntOGen archives. Writes only to
output/paper1_intogen_independence_audit/. No manuscript file is touched; nothing is committed to any Paper 1 repository.

RULE LOGIC (from intogen-plus tag v2024, commit 29e9ec5f; see intogen_filtering_rules.md):
  global q  = QVALUE_COMBINATION (weighted Stouffer over all genes)           TIER in {1,3}  <=>  q <= 0.05 (pipeline passes the .05 tier output)
  TIER 2    = CGC rescue only (q computed among CGC genes alone < 0.25) -> row can exist only for a CGC gene
  driver    = (TIER<=3 and CGC gene)  OR  (>=2 significant methods and TIER<=3)
  filters   = ... samples with >2 mutations (non-CGC only) ... germline (non-CGC only) ... lack of literature
              (n_papers==0, non-CGC, not whitelisted)
Row-level audit classes for catalog (PASS) rows:
  A  q<=0.05 and >=2 significant methods and no CGC-only filter exemption used  (would be called by the non-CGC rule)
  B  q<=0.05 but admitted through the CGC waiver of the >=2-method requirement and/or a CGC-only filter exemption
  C  q>0.05: present only through CGC rescue (TIER 2)
Pair classes (best row among cancer-matched PASS rows): A > B > C ; D = required fields missing ; NOT_A_HIT = not in catalog.
CGC membership is NOT in the released tables; 'certain CGC' is inferred only where the rules make passing impossible for a non-CGC gene.
"""
import ast
import hashlib
import json
import os
import re
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
P1 = ROOT  # this repository
OUT = ROOT / "data/derived/intogen_independence_audit"
INT_DIR = P1 / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers"
COH = P1 / "data/external/intogen/cohorts/2024-06-18_IntOGen-Cohorts/cohorts.tsv"
WL = ROOT / "data/external/intogen_plus/white_listed.txt"  # bbglab/intogen-plus tag v2024, build/datasets/postprocess/white_listed.txt
HGNC = P1 / "reference_engine_code/Human Genes.csv"
DB_DIR = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_processed"
CANCER_DB = {"Breast Cancer": "breast", "Lung Cancer": "lung", "Colorectal Cancer": "colorectal", "Prostate Cancer": "prostate", "Melanoma": "melanoma", "Bladder Cancer": "bladder",
             "Kidney Cancer": "kidney", "Pancreatic Cancer": "pancreas", "Liver Cancer": "liver", "Stomach Cancer": "stomach", "Esophageal Cancer": "esophageal", "Ovarian Cancer": "ovarian",
             "Endometrial Cancer": "endometrial", "Cervical Cancer": "cervical", "Thyroid Cancer": "thyroid", "Brain Cancer": "brain", "Leukemia": "leukemia", "Lymphoma": "lymphoma", "Myeloma": "myeloma"}
Q_THRESH = 0.05


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def extract_map(path, var):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    for n in ast.walk(tree):
        if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == var for t in n.targets):
            return ast.literal_eval(n.value)
    raise KeyError(var)


def load_tables():
    comp = pd.read_csv(INT_DIR / "Compendium_Cancer_Genes.tsv", sep="\t", encoding="utf-8-sig", low_memory=False)
    unf = pd.read_csv(INT_DIR / "Unfiltered_drivers.tsv", sep="\t", encoding="utf-8-sig", low_memory=False)
    coh = pd.read_csv(COH, sep="\t", encoding="utf-8-sig")
    wl = set(WL.read_text().split())
    comp["SYMBOL"] = comp.SYMBOL.str.upper()
    unf["SYMBOL"] = unf.SYMBOL.str.upper()
    unf["nb"] = unf.SIG_METHODS.fillna("").map(lambda s: 0 if s in ("", "combination") else len([x for x in s.split(",") if x]))
    q = unf.QVALUE_COMBINATION
    unf["q_le_005"] = q <= Q_THRESH
    unf["MUTS_PER_SAMPLE"] = unf.MUTATIONS / unf.SAMPLES_COHORT
    unf["np_missing"] = unf.NUM_PAPERS.isna()
    unf["np_zero"] = unf.NUM_PAPERS == 0
    unf["gene_whitelisted"] = unf.SYMBOL.isin(wl)
    # row-level facts that hold ONLY for CGC genes (rule logic)
    unf["cgc_waiver_methods"] = (unf.nb < 2)                                   # admitted with <2 significant methods
    unf["cgc_rescue_tier2"] = (unf.TIER == 2)
    unf["cgc_exempt_samples_filter"] = unf.SAMPLES_3MUTS.fillna(0) >= 1
    unf["cgc_exempt_germline_filter"] = unf.WARNING_GERMLINE.astype(str).eq("True") & (unf.MUTS_PER_SAMPLE > 1.0)
    unf["lit_exempt"] = unf.np_zero & ~unf.gene_whitelisted                      # PASS here => CGC (non-CGC would be 'Lack of literature evidence')
    unf["is_pass"] = unf.FILTER == "PASS"
    p = unf.is_pass
    unf["row_class"] = ""
    uses_filter_exempt = unf.cgc_exempt_samples_filter | unf.cgc_exempt_germline_filter
    unf.loc[p & unf.q_le_005 & (unf.nb >= 2) & ~uses_filter_exempt, "row_class"] = "A"
    unf.loc[p & unf.q_le_005 & ((unf.nb < 2) | uses_filter_exempt), "row_class"] = "B"
    unf.loc[p & ~unf.q_le_005, "row_class"] = "C"
    unf.loc[p & unf.q_le_005.isna(), "row_class"] = "D"
    unf["statistical_call_noncgc_rule"] = unf.q_le_005 & (unf.nb >= 2)            # statistical criterion alone, any filter status
    unf["certain_cgc_row"] = p & (unf.cgc_waiver_methods | unf.cgc_rescue_tier2 | unf.cgc_exempt_samples_filter | unf.cgc_exempt_germline_filter | unf.lit_exempt)
    return comp, unf, coh, wl


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    comp, unf, coh, wl = load_tables()
    # --- mapping: Supplementary Table S8 as typed in Paper 1's independent verifier == mapping used by the confidence script
    map_v = extract_map(P1 / "scripts/analysis/verify_headline_results.py", "CANCER_MAP")
    map_c = extract_map(P1 / "scripts/analysis/intogen_confidence_breakdown.py", "INTOGEN_CANCER_MAP")
    map_v = {k: set(v) for k, v in map_v.items()}
    map_c = {k: set(v) for k, v in map_c.items()}
    MAP_DIFF = {k: sorted(map_c[k] ^ map_v[k]) for k in map_v if map_c[k] != map_v[k]}   # known: Leukemia CML (0 rows in release 2024.09.20)
    CMAP = map_v    # published Supplementary Table S8 (typed in Paper 1's independent verifier)
    pairs = pd.read_csv(P1 / "data/canonical/precision_at_k_detail.csv")
    pairs = pairs[pairs["rank"] <= 10].copy()
    pairs["gene_norm"] = pairs.gene.str.upper()
    assert len(pairs) == 190 and pairs.cancer.nunique() == 19
    hg = pd.read_csv(HGNC, encoding="utf-8-sig").fillna("")
    approved = set(hg.symbol2.str.upper())
    alias = {}
    for _, r in hg.iterrows():
        for a in re.split(r"[,;|]\s*", str(r["Alias symbols"])):
            if a.strip():
                alias.setdefault(a.strip().upper(), set()).add(str(r.symbol2).upper())

    # --- cohort-level rows
    code_of = {c: k for k, v in CMAP.items() for c in v}
    unf["study_cancer"] = unf.CANCER_TYPE.map(code_of)
    cohort_rows, prow = [], []
    gene_all = unf.groupby("SYMBOL")
    certain_cgc_genes = set(unf[unf.certain_cgc_row].SYMBOL)
    for _, pr in pairs.iterrows():
        g, k = pr.gene_norm, pr.cancer
        gr = gene_all.get_group(g) if g in gene_all.groups else unf.iloc[0:0]
        matched = gr[gr.study_cancer == k]
        other_pass = gr[(gr.study_cancer != k) & gr.is_pass]
        for scope, df in (("cancer_matched", matched), ("other_cancer_PASS", other_pass)):
            for _, r in df.iterrows():
                cohort_rows.append(dict(cancer=k, rank=pr["rank"], gene=g, scope=scope, COHORT=r.COHORT, CANCER_TYPE=r.CANCER_TYPE, TRANSCRIPT=r.TRANSCRIPT, MUTATIONS=r.MUTATIONS,
                                        SAMPLES_COHORT=r.SAMPLES_COHORT, ALL_METHODS=r.ALL_METHODS, SIG_METHODS=r.SIG_METHODS, n_significant_methods=r.nb, QVALUE_COMBINATION=r.QVALUE_COMBINATION,
                                        RANKING=r.RANKING, TIER=r.TIER, ROLE=r.ROLE, SAMPLES_3MUTS=r.SAMPLES_3MUTS, WARNING_GERMLINE=r.WARNING_GERMLINE, NUM_PAPERS=r.NUM_PAPERS,
                                        DRIVER=r.DRIVER, FILTER=r.FILTER, is_catalog_row=r.is_pass, row_class=r.row_class, statistical_call_noncgc_rule=bool(r.statistical_call_noncgc_rule),
                                        cgc_waiver_methods=bool(r.cgc_waiver_methods) if r.is_pass else False, cgc_rescue_tier2=bool(r.cgc_rescue_tier2) if r.is_pass else False,
                                        cgc_exempt_samples_filter=bool(r.cgc_exempt_samples_filter) if r.is_pass else False, cgc_exempt_germline_filter=bool(r.cgc_exempt_germline_filter) if r.is_pass else False,
                                        literature_exemption_used=bool(r.lit_exempt) if r.is_pass else False, np_missing=bool(r.np_missing), gene_whitelisted=bool(r.gene_whitelisted)))
        mp, op = matched[matched.is_pass], other_pass
        rank = {"A": 3, "B": 2, "C": 1, "": 0}
        best = lambda df: max(df.row_class, key=lambda c: rank[c]) if len(df) else ""
        # pair-level fields
        d = dict(cancer=k, oncodigger_rank=int(pr["rank"]), gene_submitted=pr.gene, gene_norm=g, gene_is_hgnc_approved=g in approved,
                 gene_alias_of=";".join(sorted(alias.get(g, []))) if g not in approved else "",
                 global_hit_published_definition=bool(g in set(comp.SYMBOL)), exact_hit_published_definition=bool(len(mp)),
                 matched_tumour_codes=";".join(sorted(set(matched.CANCER_TYPE))), matched_cohorts_all_rows=len(matched), matched_cohorts_PASS=len(mp),
                 matched_cohort_ids_PASS=";".join(mp.COHORT), other_cancer_PASS_cohorts=len(op), other_cancer_codes_PASS=";".join(sorted(set(op.CANCER_TYPE))),
                 best_row_class_cancer_matched=best(mp) or ("NOT_A_HIT" if not len(mp) else "D"),
                 best_row_class_any_cancer=best(pd.concat([mp, op])) or "NOT_A_HIT",
                 all_matched_PASS_rows_class_A=bool(len(mp) and (mp.row_class == "A").all()),
                 min_q_matched=float(mp.QVALUE_COMBINATION.min()) if len(mp) else np.nan, max_sig_methods_matched=int(mp.nb.max()) if len(mp) else np.nan,
                 min_tier_matched=float(mp.TIER.min()) if len(mp) else np.nan, NUM_PAPERS_gene=float(gr.NUM_PAPERS.dropna().iloc[0]) if len(gr.NUM_PAPERS.dropna()) else np.nan,
                 gene_certain_cgc_lower_bound=bool(g in certain_cgc_genes), gene_whitelisted=bool(g in wl),
                 matched_candidate_rows_not_in_catalog=len(matched) - len(mp), matched_filter_reasons_not_in_catalog=";".join(sorted(set(matched[~matched.is_pass].FILTER))),
                 any_cancer_filter_reasons_for_gene=";".join(sorted(set(gr.FILTER))) if len(gr) else "")
        prow.append(d)
    P = pd.DataFrame(prow)
    C = pd.DataFrame(cohort_rows)
    # --- literature volume (audit-derived, from the 19 corpora, read-only)
    vol = {}
    for k, db in CANCER_DB.items():
        con = sqlite3.connect(f"file:{(DB_DIR / (db + '_cancer_inclusive_all') / 'lexical_index.db').as_posix()}?mode=ro", uri=True)
        mut = set(r[0] for r in con.execute("select distinct pmid from postings where token in ('mutation','mutations')"))
        for g in P[P.cancer == k].gene_norm:
            gp = set(r[0] for r in con.execute("select distinct pmid from postings where token=?", (g.lower(),)))
            vol[(k, g)] = (len(gp), len(gp & mut))
        con.close()
    P["corpus_docs_with_gene_token"] = [vol[(c, g)][0] for c, g in zip(P.cancer, P.gene_norm)]
    P["corpus_docs_gene_and_mutation_token"] = [vol[(c, g)][1] for c, g in zip(P.cancer, P.gene_norm)]
    # --- categories
    def cat(row_class, hit):
        if not hit:
            return "NOT_BENCHMARK_HIT"
        return {"A": "A", "B": "B", "C": "C", "D": "D", "": "D"}[row_class]
    P["category_cancer_specific_published_hits"] = [cat(c, h) for c, h in zip(P.best_row_class_cancer_matched, P.exact_hit_published_definition)]
    P["category_any_cancer_published_hits"] = [cat(c, h) for c, h in zip(P.best_row_class_any_cancer, P.global_hit_published_definition)]
    P.to_csv(OUT / "paper1_190_pair_evidence.csv", index=False)
    C.to_csv(OUT / "mutation_statistical_support.csv", index=False)

    # --- cancer-specific audit table
    cs = P[["cancer", "oncodigger_rank", "gene_norm", "global_hit_published_definition", "exact_hit_published_definition", "matched_tumour_codes", "matched_cohorts_all_rows", "matched_cohorts_PASS",
            "matched_cohort_ids_PASS", "other_cancer_PASS_cohorts", "other_cancer_codes_PASS", "best_row_class_cancer_matched", "best_row_class_any_cancer", "all_matched_PASS_rows_class_A",
            "matched_candidate_rows_not_in_catalog", "matched_filter_reasons_not_in_catalog", "category_cancer_specific_published_hits", "category_any_cancer_published_hits"]].copy()
    cs["support_type"] = np.where(cs.exact_hit_published_definition, "cancer_matched_catalog_support",
                                  np.where(cs.global_hit_published_definition, "global_only_gene_is_driver_in_other_cancer", "not_in_catalog"))
    cs.to_csv(OUT / "cancer_specific_support_audit.csv", index=False)

    # --- cancermine dependency audit
    D = P[["cancer", "oncodigger_rank", "gene_norm", "global_hit_published_definition", "exact_hit_published_definition", "NUM_PAPERS_gene", "gene_whitelisted", "gene_certain_cgc_lower_bound"]].copy()
    mp_rows = C[(C.scope == "cancer_matched") & C.is_catalog_row]
    gp = C[C.is_catalog_row].groupby(["cancer", "gene"])
    D["catalog_rows_any_scope"] = [len(gp.get_group((c, g))) if (c, g) in gp.groups else 0 for c, g in zip(D.cancer, D.gene_norm)]
    D["n_papers_missing_in_release"] = D.NUM_PAPERS_gene.isna()
    D["has_cancermine_papers"] = D.NUM_PAPERS_gene > 0
    D["catalog_rows_with_literature_exemption"] = [int(C[(C.cancer == c) & (C.gene == g) & C.is_catalog_row].literature_exemption_used.sum()) for c, g in zip(D.cancer, D.gene_norm)]
    flt = P.matched_filter_reasons_not_in_catalog + ";" + P.any_cancer_filter_reasons_for_gene
    D["gene_has_rows_removed_by_literature_filter"] = flt.str.contains("Lack of literature")
    D["would_become_hit_if_literature_filter_removed_exact"] = (~D.exact_hit_published_definition) & P.matched_filter_reasons_not_in_catalog.str.contains("Lack of literature")
    D["would_become_hit_if_literature_filter_removed_global"] = (~D.global_hit_published_definition) & P.any_cancer_filter_reasons_for_gene.str.contains("Lack of literature")
    D["hit_lost_if_literature_filter_removed"] = False   # by construction: the filter only removes candidates (v2024 perform_vetting)
    D["classification"] = np.select([~D.global_hit_published_definition, D.has_cancermine_papers & ~D.gene_certain_cgc_lower_bound, D.has_cancermine_papers & D.gene_certain_cgc_lower_bound,
                                     D.NUM_PAPERS_gene == 0, D.n_papers_missing_in_release],
                                    ["not_a_hit", "has_papers_literature_condition_satisfied_CGC_status_unknown", "has_papers_gene_is_certain_CGC", "zero_papers_hit_via_CGC_or_whitelist_exemption", "papers_missing_in_release"], "other")
    D.to_csv(OUT / "cancermine_dependency_audit.csv", index=False)
    print(P.category_cancer_specific_published_hits.value_counts().to_dict(), P.category_any_cancer_published_hits.value_counts().to_dict())
    print(D.classification.value_counts().to_dict())
    return P, C, D, comp, unf, CMAP


if __name__ == "__main__":
    main()
