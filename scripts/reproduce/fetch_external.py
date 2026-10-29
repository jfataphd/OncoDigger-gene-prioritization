"""
fetch_external.py -- download the CancerMine files the comparator scripts read (data/external/cancermine/) from their fixed, versioned
Zenodo record 7689627 (CancerMine, published 2023-03-01) and verify them against the SHA-256 recorded in data/external/PROVENANCE.md.
Nothing is overwritten if the existing file already matches. CancerMine is the only external comparator that has a stable public record;
the DISEASES and PubTator subsets are described in data/external/PROVENANCE.md and docs/REPRODUCTION.md.
Run: python -X utf8 scripts/reproduce/fetch_external.py
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

DEST = Path(__file__).resolve().parents[2] / "data" / "external" / "cancermine"
FILES = {  # name: (zenodo md5, sha256 recorded in PROVENANCE.md)
    "cancermine_collated.tsv": ("03a9277eb87f6f7a4d145b8104cbef06", "3f5bfc625567439cb0e2e1f96e89d96b149d2d1070457cff1c464563da233638"),
    "cancermine_sentences.tsv": ("defa8f4e07dc8e1400ba363284f562c6", "aa4fd0901d343bdf2987a76a2d60ca8a39d58d73c87333a2b750393403e6a74d"),
}


def digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    bad = 0
    for name, (md5, sha) in FILES.items():
        p = DEST / name
        if not (p.exists() and digest(p, "sha256") == sha):
            print("downloading", name, flush=True)
            urllib.request.urlretrieve(f"https://zenodo.org/api/records/7689627/files/{name}/content", p)
        ok = digest(p, "sha256") == sha and digest(p, "md5") == md5
        print(f"{name}: {'OK, matches the recorded checksums' if ok else 'CHECKSUM MISMATCH'}")
        bad += not ok
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
