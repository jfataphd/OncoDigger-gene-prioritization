"""
fetch_secondary_pmids.py -- PMIDs indexed under a collection's cancer MeSH term with the /secondary (metastasis) subheading.

Why: each collection is built with `(<MeSH term>[mh]) AND english[la] AND hasabstract`, and [mh] covers every subheading of the term, so the
organ-level collections contain papers about metastases to that organ (Liver, Brain, Lung, Kidney most), and the Melanoma term's /secondary papers concern metastatic melanoma itself. A reviewer asked for the share
of such records and for a rerun without them. This script lists them (public PubMed E-utilities, yearly windows so that no query returns more
than 10,000 records), intersects them with the frozen collection, and writes one PMID list per collection and a summary.

Run:    ONCODIGGER_ROOT=<frozen 2026-09-27 root> python -X utf8 scripts/analysis/fetch_secondary_pmids.py [--force]   (all 19 collections; existing lists are kept unless --force)
Output: data/derived/metastatic_exclusion/secondary_pmids_<collection>.csv (pmid, in_frozen_collection) and secondary_summary.csv
The PubMed counts are those of the day the script is run; records indexed after 2026-09-27 are outside the frozen build and are reported but not removed.
"""
import csv
import json
import os
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ONC = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT = ROOT / "data" / "derived" / "metastatic_exclusion"
# collection folder key -> (PubMed MeSH term used to build it, label)
COLLECTIONS = {
    "breast": "breast neoplasms", "lung": "lung neoplasms", "colorectal": "colorectal neoplasms", "prostate": "prostatic neoplasms",
    "melanoma": "melanoma", "bladder": "urinary bladder neoplasms", "kidney": "kidney neoplasms", "pancreas": "pancreatic neoplasms",
    "liver": "liver neoplasms", "stomach": "stomach neoplasms", "esophageal": "esophageal neoplasms", "ovarian": "ovarian neoplasms",
    "endometrial": "endometrial neoplasms", "cervical": "uterine cervical neoplasms", "thyroid": "thyroid neoplasms",
    "brain": "brain neoplasms", "leukemia": "leukemia", "lymphoma": "lymphoma", "myeloma": "multiple myeloma",
}
BASE = " AND english[la] AND hasabstract"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def esearch(term: str, mindate: str, maxdate: str) -> list[str]:
    params = {"db": "pubmed", "retmode": "json", "retmax": "10000", "datetype": "pdat", "mindate": mindate, "maxdate": maxdate,
              "term": term}
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            res = json.load(urllib.request.urlopen(url, timeout=90))["esearchresult"]
            ids, count = res["idlist"], int(res["count"])
            if count > len(ids):
                raise RuntimeError(f"window {mindate}-{maxdate} has {count} records, more than the 10,000 returned")
            return ids
        except RuntimeError:
            raise
        except Exception:
            time.sleep(3 * (attempt + 1))
    raise RuntimeError("PubMed request failed: " + url)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    force = "--force" in sys.argv
    for key, mesh in COLLECTIONS.items():
        if (OUT / f"secondary_pmids_{key}.csv").exists() and not force:
            continue
        term = f'("{mesh}/secondary"[mh])' + BASE
        pmids = []
        for year in range(2000, 2027):
            hi = f"{year}/12/31" if year < 2026 else "2026/09/27"
            pmids += esearch(term, f"{year}/01/01", hi)
            time.sleep(0.4)
        pmids = sorted(set(pmids), key=int)
        con = sqlite3.connect(f"file:{(ONC / 'data' / 'pubmed_processed' / f'{key}_cancer_inclusive_all' / 'lexical_index.db').as_posix()}?mode=ro&immutable=1", uri=True)
        frozen = {str(r[0]) for r in con.execute("select pmid from documents")}
        con.close()
        with (OUT / f"secondary_pmids_{key}.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["pmid", "in_frozen_collection"])
            for p in pmids:
                w.writerow([p, p in frozen])
        inside = sum(p in frozen for p in pmids)
        summary.append({"collection": key, "mesh_term": mesh, "pubmed_secondary_records": len(pmids), "in_frozen_collection": inside,
                        "frozen_collection_size": len(frozen), "share_of_frozen_collection_pct": round(100 * inside / len(frozen), 2)})
        print(summary[-1], flush=True)
    prev = list(csv.DictReader((OUT / "secondary_summary.csv").open(encoding="utf-8"))) if (OUT / "secondary_summary.csv").exists() else []
    allrows = [r for r in prev if r["collection"] not in {x["collection"] for x in summary}] + summary
    with (OUT / "secondary_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(allrows[0]))
        w.writeheader()
        w.writerows(allrows)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
