"""
rebuild_corpus_from_pmids.py -- rebuild one collection's lexical_index.db from NCBI, starting from the exact PMID list the paper used.

The paper's 19 collections are not redistributed (they hold PubMed abstracts). Each collection's PMID list is in
data/provenance/pmid_manifests/<collection>.pmids.gz. This script downloads exactly those records from NCBI E-utilities
(efetch, by PMID, so the paper set is fixed even though PubMed has grown), cleans and indexes them with the same steps as the
application's pipeline (scripts/pubmed_fetch.py parse, scripts/pubmed_preprocess.py clean, build_lexical_index.py +
migrate_to_sqlite.py schema; reproduced here), and writes lexical_index.db in the same schema the ranking scripts read.

PubMed records are revised after the fact (abstract corrections, late MeSH indexing, withdrawals), so a rebuilt database is not
byte-identical to the author's. --compare points at a stored database and reports exactly how many records differ in title,
abstract, MeSH, year or length, so the drift is measured and not assumed.

Usage:
  python -X utf8 scripts/reproduce/rebuild_corpus_from_pmids.py myeloma --out <dir> [--limit N] [--email you@x.org] [--api-key KEY]
                                                              [--compare <stored lexical_index.db>]
Needs network access to eutils.ncbi.nlm.nih.gov. Without an API key NCBI allows 3 requests/second (about 200 records each);
the script paces itself. Raw downloads are cached in <out>/raw_<collection>.jsonl so an interrupted run resumes.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import math
import re
import sqlite3
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
LEX = ROOT / "reference_engine_code"
TOKEN_RE = re.compile(r"[a-z0-9]+(?:[-'][a-z0-9]+)*")
SPACE_RE = re.compile(r"\s+")

# ---- fetch and parse (as scripts/pubmed_fetch.py) -------------------------------------------------------------------------


def text_of(el):
    return "" if el is None else " ".join(p.strip() for p in el.itertext() if p and p.strip())


def pub_date(article):
    el = article.find(".//JournalIssue/PubDate")
    if el is None:
        return ""
    y, m, d, med = (text_of(el.find(k)) for k in ("Year", "Month", "Day", "MedlineDate"))
    return "-".join(p for p in (y, m, d) if p) or med


def parse_articles(xml_text: str) -> list[dict]:
    out = []
    for cit in ET.fromstring(xml_text).findall(".//PubmedArticle"):
        art, med = cit.find(".//Article"), cit.find(".//MedlineCitation")
        if art is None or med is None:
            continue
        ids = {}
        for a in cit.findall(".//PubmedData/ArticleIdList/ArticleId"):
            if a.attrib.get("IdType") and text_of(a):
                ids[a.attrib["IdType"]] = text_of(a)
        authors = []
        for a in art.findall(".//AuthorList/Author"):
            authors.append({"last_name": text_of(a.find("LastName")), "fore_name": text_of(a.find("ForeName")),
                            "initials": text_of(a.find("Initials")), "collective_name": text_of(a.find("CollectiveName")),
                            "affiliations": []})
        out.append({
            "pmid": text_of(med.find("PMID")), "title": text_of(art.find("ArticleTitle")),
            "abstract": " ".join(text_of(p) for p in art.findall(".//Abstract/AbstractText")).strip(),
            "journal": text_of(art.find(".//Journal/Title")), "journal_iso": text_of(art.find(".//Journal/ISOAbbreviation")),
            "publication_date": pub_date(art), "article_ids": ids, "authors": authors,
            "publication_types": [text_of(p) for p in art.findall(".//PublicationType")],
            "mesh_terms": [text_of(m.find("DescriptorName")) for m in med.findall(".//MeshHeading")],
        })
    return out


def efetch(pmids: list[str], email, api_key, retries=5) -> str:
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(EUTILS, data=data), timeout=120) as r:
                return r.read().decode("utf-8")
        except Exception as e:  # noqa: BLE001
            wait = 2 * (attempt + 1)
            print(f"  fetch error ({e}); retry in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise SystemExit("NCBI fetch failed repeatedly; rerun to resume from the cached records")


def download(pmids, raw_path: Path, email, api_key, batch=200):
    have = set()
    if raw_path.exists():
        with raw_path.open(encoding="utf-8") as f:
            have = {json.loads(line)["pmid"] for line in f if line.strip()}
    todo = [p for p in pmids if p not in have]
    pause = 0.11 if api_key else 0.35
    print(f"{len(have):,} cached, {len(todo):,} to fetch in {math.ceil(len(todo) / batch):,} requests")
    with raw_path.open("a", encoding="utf-8") as out:
        for i in range(0, len(todo), batch):
            for rec in parse_articles(efetch(todo[i:i + batch], email, api_key)):
                out.write(json.dumps(rec) + "\n")
            out.flush()
            if (i // batch) % 20 == 0:
                print(f"  {min(i + batch, len(todo)):,}/{len(todo):,}", flush=True)
            time.sleep(pause)


# ---- clean (as scripts/pubmed_preprocess.py) --------------------------------------------------------------------------------


def normalize_text(text: str) -> str:
    text = html.unescape(text or "")
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii").lower()
    text = text.replace("/", " ").replace("_", " ")
    text = re.sub(r"[‐-―]", "-", text)
    return SPACE_RE.sub(" ", re.sub(r"[^a-z0-9\-\s']", " ", text)).strip()


def csv_column(path: Path, column: str) -> list[str]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return [row[column].strip() for row in csv.DictReader(f) if row.get(column, "").strip()]


def phrase_replacer():
    terms = list(csv_column(LEX / "Human Genes.csv", "Approved name"))
    terms += csv_column(LEX / "kegg_drug_list_lowercase.csv", "drugs")
    comp = csv_column(LEX / "kegg_compounds_lowercase.csv", "compound")
    for value in comp:  # split_compound_synonyms
        terms += [p.strip() for p in value.split(";") if p.strip()]
    terms += csv_column(LEX / "phytochemicals.csv", "phyto")
    norm = {normalize_text(t) for t in terms}
    phrases = sorted((t for t in norm if " " in t and len(t) >= 5), key=len, reverse=True)
    trie: dict = {}
    for ph in phrases:
        node = trie
        for ch in ph:
            node = node.setdefault(ch, {})
        node[None] = {}

    def expr(node):
        branches = [re.escape(c) + expr(ch) for c, ch in sorted((c, n) for c, n in node.items() if c is not None)]
        if None in node:
            branches.append("")
        return branches[0] if len(branches) == 1 else "(?:" + "|".join(branches) + ")"

    pat = re.compile(r"(?<![a-z0-9-])(" + expr(trie) + r")(?![a-z0-9-])")
    return lambda text: pat.sub(lambda m: m.group(1).replace(" ", "-"), text)


def year_of(date: str):
    m = re.search(r"\b(19|20)\d{2}\b", date or "")
    return int(m.group(0)) if m else None


def clean(records, min_tokens=20):
    replace, seen, kept = phrase_replacer(), set(), []
    for r in records:
        if r["pmid"] in seen:
            continue
        seen.add(r["pmid"])
        tokens = TOKEN_RE.findall(replace(normalize_text(r.get("abstract", ""))))
        if len(tokens) < min_tokens:
            continue
        kept.append({
            "pmid": r["pmid"], "title": normalize_text(r.get("title", "")), "abstract": normalize_text(r.get("abstract", "")),
            "publication_date": r.get("publication_date", ""), "publication_year": year_of(r.get("publication_date", "")),
            "journal": normalize_text(r.get("journal", "")), "journal_iso": normalize_text(r.get("journal_iso", "")),
            "article_ids": r.get("article_ids", {}),
            "authors": [{"last_name": normalize_text(a.get("last_name", "")), "fore_name": normalize_text(a.get("fore_name", "")),
                         "collective_name": normalize_text(a.get("collective_name", ""))} for a in r.get("authors", [])[:20]],
            "publication_types": r.get("publication_types", []),
            "mesh_terms": [normalize_text(t) for t in r.get("mesh_terms", [])], "tokens": tokens})
    return kept


# ---- index (as scripts/build_lexical_index.py + migrate_to_sqlite.py) ----------------------------------------------------------

SCHEMA = """
CREATE TABLE corpus_stats (avg_doc_length REAL NOT NULL, doc_count INTEGER NOT NULL);
CREATE TABLE idf (token TEXT PRIMARY KEY, idf REAL NOT NULL);
CREATE TABLE documents (pmid TEXT PRIMARY KEY, title TEXT, abstract TEXT, journal TEXT, journal_iso TEXT, publication_year INTEGER,
  doc_length INTEGER, mesh_terms TEXT, publication_types TEXT, authors TEXT, article_ids TEXT, publication_date TEXT);
CREATE TABLE postings (token TEXT NOT NULL, pmid TEXT NOT NULL, tf INTEGER NOT NULL);
"""


def build_db(docs, db_path: Path):
    postings, total = defaultdict(list), 0
    for d in docs:
        total += len(d["tokens"])
        for tok, n in Counter(d["tokens"]).items():
            postings[tok].append((d["pmid"], n))
    n_docs = len(docs)
    idf = {t: math.log(1 + (n_docs - len(p) + 0.5) / (len(p) + 0.5)) for t, p in postings.items()}
    if db_path.exists():
        db_path.unlink()
    con = sqlite3.connect(db_path)
    con.executescript(SCHEMA)
    con.execute("INSERT INTO corpus_stats VALUES (?,?)", (total / n_docs if n_docs else 0, n_docs))
    con.executemany("INSERT INTO idf VALUES (?,?)", idf.items())
    con.executemany("INSERT INTO documents VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", [
        (d["pmid"], d["title"], d["abstract"], d["journal"], d["journal_iso"], d["publication_year"], len(d["tokens"]),
         json.dumps(d["mesh_terms"]), json.dumps(d["publication_types"]), json.dumps(d["authors"]), json.dumps(d["article_ids"]),
         d["publication_date"]) for d in docs])
    con.executemany("INSERT INTO postings VALUES (?,?,?)", ((t, p, n) for t, rows in postings.items() for p, n in rows))
    con.execute("CREATE INDEX idx_postings_token ON postings(token)")
    con.commit()
    con.close()


def compare(new_db: Path, stored_db: Path, pmids: set[str]):
    def load(p):
        con = sqlite3.connect(f"file:{Path(p).resolve().as_posix()}?mode=ro&immutable=1", uri=True)
        rows = {r[0]: r[1:] for r in con.execute(
            "SELECT pmid, title, abstract, mesh_terms, publication_year, doc_length FROM documents")}
        con.close()
        return rows
    a, b = load(new_db), load(stored_db)
    b = {p: v for p, v in b.items() if p in pmids}  # restrict to the PMIDs asked for (matters with --limit)
    both = a.keys() & b.keys()
    fields = ["title", "abstract", "mesh_terms", "publication_year", "doc_length"]
    diff = {f: sum(1 for p in both if a[p][i] != b[p][i]) for i, f in enumerate(fields)}
    any_diff = sum(1 for p in both if a[p] != b[p])
    print(f"\nComparison with the stored database {stored_db.name}:")
    print(f"  PMIDs in the paper's list:        {len(pmids):,}")
    print(f"  in stored database:               {len(b):,}")
    print(f"  in rebuilt database:              {len(a):,}")
    print(f"  missing from the rebuild:         {len(b.keys() - a.keys()):,}  (withdrawn, now without abstract, or under min tokens)")
    print(f"  only in the rebuild:              {len(a.keys() - b.keys()):,}")
    print(f"  records present in both:          {len(both):,}")
    print(f"  of those, identical in all five fields: {len(both) - any_diff:,} ({100 * (len(both) - any_diff) / max(len(both), 1):.2f}%)")
    for f, n in diff.items():
        print(f"    differing {f:<17}{n:,}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("collection")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--email")
    ap.add_argument("--api-key")
    ap.add_argument("--compare", type=Path)
    a = ap.parse_args()
    with gzip.open(ROOT / "data" / "provenance" / "pmid_manifests" / f"{a.collection}.pmids.gz", "rt") as f:
        pmids = [x.strip() for x in f if x.strip()]
    if a.limit:
        pmids = pmids[:a.limit]
    a.out.mkdir(parents=True, exist_ok=True)
    raw = a.out / f"raw_{a.collection}.jsonl"
    t0 = time.time()
    download(pmids, raw, a.email, a.api_key)
    with raw.open(encoding="utf-8") as f:
        records = [json.loads(line) for line in f if line.strip()]
    wanted = set(pmids)
    docs = clean([r for r in records if r["pmid"] in wanted])
    db = a.out / f"{a.collection}_lexical_index.db"
    build_db(docs, db)
    print(f"\nWrote {db} with {len(docs):,} documents from {len(pmids):,} PMIDs in {time.time() - t0:.0f}s")
    if a.compare:
        compare(db, a.compare, wanted)


if __name__ == "__main__":
    main()
