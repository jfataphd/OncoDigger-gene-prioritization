"""
test_bg_matcher.py -- hand-built corpus with expected background counts worked out by hand from the document text (not from the matcher).
Covers: aliases that overlap (HER2 / HER-2 / ERBB2 in one paper), hyphenated and phrase-only symbols (NKX2-1, HLA-A), an alias the lexicon does not
contain (c-erbB-2), longest-match masking (EGF inside epidermal growth factor receptor), a suppressed ambiguous bare symbol (ELN), an ambiguous
common word (was), a symbol and a short alias that the vocabulary does not carry (TP53 / p53), and a PMID listed twice.
Run: PYTHONHASHSEED=0 python -X utf8 test_bg_matcher.py <ONCODIGGER_ROOT>
"""
import os, sys
from pathlib import Path
assert os.environ.get("PYTHONHASHSEED") == "0"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bg_matcher as bg
mt = bg.Matcher(Path(sys.argv[1]))
docs = [
    ("1", "HER2 amplification in breast tumours", ""),
    ("2", "HER-2 positive disease", "Treatment of HER-2 positive tumours."),
    ("3", "ERBB2 and HER2 co-expression", ""),
    ("4", "c-erbB-2 overexpression", ""),
    ("1", "HER2 amplification in breast tumours", ""),            # duplicate PMID 1: must not be counted twice
    ("6", "NKX2-1 and TTF-1 staining", ""),
    ("7", "Thyroid transcription factor TTF-1 in lung adenocarcinoma", ""),
    ("8", "HLA-A*02:01 restricted peptides", ""),
    ("9", "Epidermal growth factor receptor inhibitors", ""),
    ("10", "Epidermal growth factor signalling", ""),
    ("11", "Elastin fibres and the ELN gene", ""),
    ("12", "The ELN trial of early nephritis", ""),
    ("13", "p53 mutation in tumours", ""),
    ("14", "TP53 and p53 status", ""),
    ("15", "This was observed in all patients", "No gene is named here."),
]
b, s, n = mt.count(docs, pool_ids={"1", "2", "9"})
by = {}
for (eid, lab), v in b.items():
    by.setdefault(lab, 0); by[lab] += v
expected = {"ERBB2": 3, "NKX2-1": 2, "HLA-A": 1, "EGFR": 1, "EGF": 1, "ELN": 1, "TP53": 1}  # TP53: only doc 14. Bare "p53" (doc 13) is NOT recognised, because the vocabulary drops aliases shorter than 4 characters (lexicons.useful_gene_term). Documented limitation, not changed here.
fail = 0
for lab, want in expected.items():
    got = by.get(lab, 0)
    print(("ok   " if got == want else "FAIL "), lab, "expected", want, "got", got)
    fail += got != want
for lab in ("WAS",):
    got = by.get(lab, 0); print(("ok   " if got == 0 else "FAIL "), lab, "expected 0 got", got); fail += got != 0
print("papers counted:", n, "(expected 14, PMID 1 once)"); fail += n != 14
sp = {lab: v for (eid, lab), v in s.items()}
print("pool S (pool = PMIDs 1, 2, 9):", {k: sp.get(k, 0) for k in ("ERBB2", "EGFR")}, "(expected ERBB2 2, EGFR 1)")
fail += sp.get("ERBB2", 0) != 2 or sp.get("EGFR", 0) != 1
print("RESULT:", "PASS" if not fail else f"{fail} FAILURE(S)")
sys.exit(1 if fail else 0)
