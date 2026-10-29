"""
vocab_correction_v5_table.py -- the frozen, reviewable decision table of the v5 vocabulary correction (vocab_policy.py; Supplementary Methods S1).

Scope: every candidate of vocab_inventory_screen.py (whole-vocabulary screen of single-word aliases and ordinary-word symbols) whose gene appears in
any reported ranking or control (602 genes: main top-50 lists, pooled control, cutoff lists, both cleaned cells of Table 2) -- 119 of 165 candidates.
Replacement entries that enter a list after the correction are checked in the rebuild report.

Decision rule (uniform, independent of any benchmark; written before the rebuild):
  SUPPRESS  the term refers to the gene in <= 10% of sampled papers in every collection examined (seeded, random.seed(42); 6 papers read for
            every term, 100 screened with a gene-context pattern for the aliases, 14-30 read for the symbols and mixed terms).
  RETAIN    the term is a recognised biological name of the gene (protein name, marker, abbreviation of the gene product) -- including
            names that are also dictionary words (insulin, albumin, survivin, snail, trail, nestin) -- and the sample is mostly genuine.
  RETAIN-CD the term is a real gene name that is also a common word, and its sample is mixed: genuine in some collections, other uses in
            others. A bare symbol cannot be suppressed per collection without losing real recognition, so it is kept and disclosed.
Suppression is term-level: only the stated token is removed from the stated gene; the gene keeps its approved symbol (unless listed),
its full name and its other aliases. Readings were made by a language-model agent that had seen the manuscript (not independent validation).

Output: data/derived/vocab_rescoring/vocab_policy_v5_table.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

OUT = Path(__file__).resolve().parents[2] / "data" / "derived" / "vocab_rescoring" / "vocab_policy_v5_table.csv"

# kind, term, gene, decision, collection read, genuine/read (manual), meaning of the non-gene uses / reason, main top-50 positions (corrected lists)
A, S = "alias", "symbol"
SUP, RET, RCD = "suppress", "retain", "retain-collection-dependent"
ROWS = [
    # --- v4 (earlier in this correction) ---
    (S, "SHE", "SHE", SUP, "breast/thyroid", "0/8", "the pronoun", "Thyroid 49"),
    (S, "RHO", "RHO", SUP, "stomach", "0/30", "Rho GTPases, Spearman's rho", "Stomach 46"),
    (S, "CCS", "CCS", SUP, "colorectal", "3/30", "Cronkhite-Canada syndrome, colon cancers, screening and scoring systems", "Colorectal 50"),
    (S, "NODAL", "NODAL", SUP, "thyroid, lymphoma", "0/60", "the adjective nodal", "Thyroid 21; Lymphoma 24"),
    (A, "rank", "TNFRSF11A", SUP, "breast", "0/14", "log-rank, rank-sum, Spearman rank", "Lung 18; Colorectal 18; Breast 22; Prostate 40; Melanoma 48"),
    (A, "nude", "NDE1", SUP, "lung", "0/30", "nude mice", "-"),
    (A, "lobe", "AKT1S1", SUP, "lung, liver", "0/30", "anatomical lobe", "-"),
    # --- v5 aliases: suppress ---
    (A, "blast", "CD48", SUP, "leukemia", "0/6", "blast cells, blast crisis", "Leukemia 23"),
    (A, "race", "AMACR", SUP, "breast", "0/6", "race/ethnicity", "Endometrial 39; Breast 41; Kidney 45; Prostate 48"),
    (A, "flap", "ALOX5AP", SUP, "breast", "0/6", "surgical flap (reconstruction)", "-"),
    (A, "icis", "MTUS1", SUP, "lung", "0/6", "immune checkpoint inhibitors (ICIs)", "Lung 16; Bladder 50"),
    (A, "icis", "MTUS2", SUP, "lung", "0/6", "immune checkpoint inhibitors (ICIs)", "-"),
    (A, "killer", "TNFRSF10B", SUP, "lymphoma", "0/6", "natural killer cells", "-"),
    (A, "huh7", "MIR7-3HG", SUP, "liver", "0/6", "the Huh7 hepatoma cell line", "-"),
    (A, "pilot", "EGR3", SUP, "breast", "0/6", "pilot study", "Breast 45"),
    (A, "ca125", "NBR1", SUP, "ovarian", "0/6", "the CA125 tumour marker (MUC16 keeps the alias)", "Ovarian 36"),
    (A, "fact", "SUPT16H", SUP, "breast", "0/6", "'the fact that', FACT (Functional Assessment of Cancer Therapy)", "-"),
    (A, "minor", "NR4A3", SUP, "breast", "0/6", "minor complications, minor allele", "Lung 23; Lymphoma 47"),
    (A, "task", "KCNK3", SUP, "breast", "0/6", "task (cognitive, classification)", "-"),
    (A, "ctcl", "TSPYL2", SUP, "lymphoma", "0/6", "cutaneous T-cell lymphoma (CTCL)", "-"),
    (A, "maps", "C3ORF62", SUP, "brain", "0/6", "imaging parameter maps", "-"),
    (A, "larc", "CCL20", SUP, "colorectal", "0/6", "locally advanced rectal cancer (LARC)", "-"),
    (A, "ct26", "DDX53", SUP, "colorectal", "0/6", "the CT26 colon carcinoma mouse cell line", "-"),
    (A, "fast", "FASTK", SUP, "breast", "0/6", "the adjective fast, FAST", "Colorectal 28"),
    (A, "fast", "FOXD3-AS1", SUP, "breast", "0/6", "the adjective fast, FAST", "Colorectal 29; Melanoma 40; Lung 41"),
    (A, "particle", "PARTICL", SUP, "breast", "0/6", "particle size (nanoparticles)", "-"),
    (A, "nets", "SPINK5", SUP, "pancreas", "0/6", "neuroendocrine tumors (NETs)", "-"),
    (A, "cava", "CA5A", SUP, "kidney", "0/6", "inferior vena cava", "-"),
    (A, "anova", "NOVA2", SUP, "breast", "0/6", "the ANOVA test", "-"),
    (A, "nact", "SLC13A5", SUP, "breast", "0/6", "neoadjuvant chemotherapy (NACT)", "-"),
    (A, "scar", "RPS4X", SUP, "breast", "0/6", "surgical scar", "-"),
    (A, "tramp", "TNFRSF25", SUP, "prostate", "0/6", "the TRAMP transgenic mouse model", "-"),
    (A, "tramp", "TRAM1", SUP, "prostate", "0/6", "the TRAMP transgenic mouse model", "-"),
    (A, "tramp", "ZBED1", SUP, "prostate", "0/6", "the TRAMP transgenic mouse model", "-"),
    (A, "tram", "TRAM1", SUP, "breast", "0/6", "TRAM flap reconstruction", "-"),
    (A, "asparaginase", "ASPG", SUP, "leukemia", "0/6", "the drug L-asparaginase", "-"),
    (A, "sera", "PHGDH", SUP, "liver", "0/6", "plural of serum", "-"),
    (A, "wave", "WASF1", SUP, "breast", "1/6 (1/100 pattern)", "shear wave, E-wave, wave length", "-"),
    (A, "warts", "LATS1", SUP, "cervical", "0/6", "HPV genital warts (LATS1 keeps its symbol and full name)", "-"),
    (A, "call", "CHL1", SUP, "breast", "0/6", "the verb call, call rate", "-"),
    (A, "ralp", "SHC4", SUP, "prostate", "0/6", "robot-assisted laparoscopic prostatectomy (RALP)", "-"),
    (A, "rarc", "RARG", SUP, "bladder", "0/6", "robot-assisted radical cystectomy (RARC)", "-"),
    (A, "fell", "STAB2", SUP, "breast", "0/6", "the verb fell", "-"),
    (A, "spin", "GTF2I", SUP, "brain", "0/6", "MRI spin echo, arterial spin labeling", "-"),
    (A, "cats", "PIMREG", SUP, "lymphoma", "0/6", "the animal (feline lymphoma)", "-"),
    (A, "nail", "CD244", SUP, "melanoma", "0/6", "nail (subungual melanoma)", "-"),
    (A, "clip", "POMC", SUP, "breast", "0/6", "surgical clip, clip art", "Brain 47"),
    (A, "clip", "CD74", SUP, "breast", "0/6", "surgical clip, clip art", "-"),
    (A, "delta", "YY1", SUP, "leukemia", "0/6", "T-cell receptor delta, PKC delta, Delta notation", "Pancreatic 40; Cervical 41"),
    (A, "crest", "SS18L1", SUP, "melanoma", "0/6", "neural crest", "-"),
    (A, "genesis", "FOXD3", SUP, "breast", "0/6", "tumor genesis, genesis of cancer", "Cervical 30; Thyroid 43; Esophageal 47"),
    (A, "wish", "NCKIPSD", SUP, "breast", "0/6", "the verb wish", "-"),
    (A, "hits", "FAM107B", SUP, "breast", "0/6", "screen hits, second hits", "Colorectal 47"),
    (A, "cd25", "ISG20", SUP, "leukemia", "0/6", "CD25 marker (IL2RA keeps the alias)", "-"),
    (A, "nrf2", "GABPA", SUP, "lung", "0/6", "NRF2 = NFE2L2 (NFE2L2 keeps the alias)", "Esophageal 39"),
    # --- v5 symbols: suppress ---
    (S, "ESD", "ESD", SUP, "stomach, colorectal", "0/6 (0/180 pattern)", "endoscopic submucosal dissection", "-"),
    (S, "GEM", "GEM", SUP, "pancreas", "0/6 (1/180)", "gemcitabine (GEM)", "-"),
    (S, "CAD", "CAD", SUP, "breast", "0/6 (2/180)", "computer-aided diagnosis, coronary artery disease", "-"),
    (S, "ARC", "ARC", SUP, "prostate, lung, brain", "0/6 (0/180)", "volumetric modulated arc therapy", "-"),
    (S, "NHS", "NHS", SUP, "breast", "0/6 (1/180)", "National Health Service, Nurses' Health Study", "-"),
    (S, "COIL", "COIL", SUP, "prostate", "0/6 (1/180)", "endorectal and phased-array MRI coils", "-"),
    (S, "PTS", "PTS", SUP, "breast", "0/6 (0/180)", "patients (pts), phyllodes tumors", "-"),
    (S, "POR", "POR", SUP, "colorectal", "0/6 (0/60)", "Spanish 'por', postoperative recurrence, odds ratio", "-"),
    (S, "MOS", "MOS", SUP, "lung", "0/6 (2/180)", "months, median overall survival (mOS)", "-"),
    (S, "ADM", "ADM", SUP, "breast", "0/6 (5/180)", "acellular dermal matrix", "-"),
    (S, "IDE", "IDE", SUP, "liver", "0/6 (1/60)", "nucleos(t)ide analogs (tokenization)", "-"),
    (S, "IDS", "IDS", SUP, "ovarian", "0/6 (1/114)", "interval debulking surgery", "-"),
    (S, "MAX", "MAX", SUP, "lung, brain", "1/27", "SUVmax, Vmax, Cmax, Dmax", "-"),
    (S, "REST", "REST", SUP, "breast, brain", "1/28", "'the rest of', rest (at rest)", "-"),
    # --- retained: genuine names ---
    (A, "her2", "ERBB2", RET, "breast", "6/6", "-", "-"), (A, "aromatase", "CYP19A1", RET, "breast", "6/6", "aromatase inhibitors (the enzyme)", "-"),
    (A, "psma", "FOLH1", RET, "prostate", "6/6", "-", "-"), (A, "albumin", "ALB", RET, "liver", "6/6", "protein name; serum albumin", "-"),
    (A, "parp", "PARP1", RET, "breast", "6/6", "family name (PARP1/2)", "-"), (A, "tyrosinase", "TYR", RET, "melanoma", "6/6", "-", "-"),
    (A, "insulin", "INS", RET, "pancreas", "6/6", "protein name", "-"), (A, "pi3k", "PIK3CA", RET, "liver", "6/6", "family name", "-"),
    (A, "ca125", "MUC16", RET, "ovarian", "6/6", "-", "-"), (A, "vimentin", "VIM", RET, "lung", "6/6", "-", "-"),
    (A, "cd10", "MME", RET, "lymphoma", "6/6", "-", "-"), (A, "somatostatin", "SST", RET, "pancreas", "6/6", "-", "-"),
    (A, "erk1", "MAPK3", RET, "breast", "6/6", "-", "-"), (A, "cd56", "NCAM1", RET, "lymphoma", "6/6", "-", "-"),
    (A, "acth", "POMC", RET, "brain", "6/6", "hormone product of the gene", "Brain 47"), (A, "cd138", "SDC1", RET, "myeloma", "6/6", "-", "-"),
    (A, "survivin", "BIRC5", RET, "breast", "6/6", "-", "-"), (A, "cd133", "PROM1", RET, "colorectal", "6/6", "-", "-"),
    (A, "snail", "SNAI1", RET, "breast", "6/6", "the transcription factor", "-"), (A, "gp100", "PMEL", RET, "melanoma", "6/6", "-", "-"),
    (A, "cd25", "IL2RA", RET, "leukemia", "6/6", "-", "-"), (A, "leptin", "LEP", RET, "breast", "6/6", "-", "-"),
    (A, "nrf2", "NFE2L2", RET, "lung", "6/6", "-", "-"), (A, "hdac", "HDAC9", RET, "leukemia", "6/6", "family name (not specific to HDAC9; disclosed)", "-"),
    (A, "hsp90", "HSP90AA1", RET, "breast", "4/4", "family name", "-"), (A, "trail", "TNFSF10", RET, "colorectal", "6/6", "-", "-"),
    (A, "nestin", "NES", RET, "brain", "6/6", "-", "-"), (A, "gastrin", "GAST", RET, "stomach", "6/6", "-", "-"),
    (A, "synaptophysin", "SYP", RET, "brain", "6/6", "-", "-"), (A, "cd23", "FCER2", RET, "lymphoma", "6/6", "-", "-"),
    (A, "ck20", "KRT20", RET, "colorectal", "6/6", "-", "-"), (A, "ck19", "KRT19", RET, "liver", "6/6", "-", "-"),
    (A, "rankl", "TNFSF11", RET, "breast", "6/6", "-", "-"), (A, "vegfr", "KDR", RET, "kidney", "6/6", "family name (VEGFR1-3; disclosed)", "-"),
    (A, "adiponectin", "ADIPOQ", RET, "breast", "6/6", "-", "-"), (A, "glucagon", "GCG", RET, "pancreas", "6/6", "-", "-"),
    (A, "hmb45", "PMEL", RET, "melanoma", "6/6", "antibody to gp100", "-"), (A, "cd43", "SPN", RET, "lymphoma", "6/6", "-", "-"),
    (A, "ptc1", "CCDC6", RET, "thyroid", "6/6", "RET/PTC1 fusion partner", "-"),
    (S, "AFP", "AFP", RET, "liver", "6/6", "abbreviation of alpha-fetoprotein", "Liver 38"), (S, "MGMT", "MGMT", RET, "brain", "6/6", "-", "Brain 9; Esophageal 49"),
    (S, "RET", "RET", RET, "thyroid, lung", "6/6", "-", "Thyroid 3; Lung 20; Cervical 45; Colorectal 49"),
    (S, "ATM", "ATM", RET, "breast", "6/6", "-", "16 lists"), (S, "TYR", "TYR", RET, "melanoma", "6/6", "-", "Thyroid 33"),
    (S, "SRC", "SRC", RET, "breast, leukemia", "5/6", "-", "-"),
    # --- retained: real gene names that are also common words, collection-dependent ---
    (A, "chop", "DDIT3", RCD, "breast/lung/pancreas/liver; lymphoma", "8/14, 11/14, 11/14, 13/14; 0/6", "CHOP regimen in lymphoma and case reports; ER-stress gene elsewhere", "Lymphoma 50"),
    (S, "KIT", "KIT", RCD, "colorectal; leukemia", "1/14; 8/14", "assay kit, FIT/screening kit", "14 lists (Melanoma 3 ... Bladder 42)"),
    (S, "MET", "MET", RCD, "breast; lung", "1/14; 6/14", "'met the criteria', metabolic equivalent (MET)", "13 lists (Lung 9 ... Endometrial 45)"),
    (S, "APC", "APC", RCD, "lung; breast; colorectal", "4/13; 7/14; 5/6", "annual percent change, age-period-cohort, anaphase-promoting complex", "11 lists (Stomach 7 ... Ovarian 41)"),
    (S, "POLE", "POLE", RCD, "kidney; endometrial", "0/14; 56/60 (pattern)", "upper/lower pole of the kidney", "Endometrial 3; Ovarian 5; Colorectal 12; Brain 18; Stomach 49; Pancreatic 50"),
    (S, "MIA", "MIA", RCD, "pancreas; lung; melanoma", "0/14; 1/14; 59/60 (pattern)", "MIA PaCa-2 cell line; minimally invasive adenocarcinoma", "-"),
    (S, "BAD", "BAD", RCD, "breast; leukemia", "3/14; 6/14", "'bad prognosis'", "-"),
    (S, "FAS", "FAS", RCD, "breast; colorectal; leukemia", "mostly genuine (33-46/60 pattern)", "fibroadenomas (FAs), fatty acid synthase", "Lymphoma 36"),
]


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "term", "gene", "decision", "collections_read", "genuine_over_read", "non_gene_meaning", "main_top50_positions"])
        w.writerows(ROWS)
    n = {d: sum(r[3] == d for r in ROWS) for d in (SUP, RET, RCD)}
    print(f"wrote {OUT}: {len(ROWS)} rows; suppress {n[SUP]} ({sum(r[3] == SUP and r[0] == A for r in ROWS)} aliases, "
          f"{sum(r[3] == SUP and r[0] == S for r in ROWS)} symbols), retain {n[RET]}, retain-collection-dependent {n[RCD]}")


if __name__ == "__main__":
    main()
