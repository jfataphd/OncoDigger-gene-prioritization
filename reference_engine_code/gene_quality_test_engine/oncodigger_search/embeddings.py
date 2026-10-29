from pathlib import Path
import json
import re

from .lexical import metadata_multiplier, snippet


TOKEN_RE = re.compile(r"[a-z0-9]+(?:[-'][a-z0-9]+)*")


def embedding_artifact_exists(index_dir: str | Path) -> bool:
    if not index_dir:
        return False
    path = Path(index_dir)
    return (path / "manifest.json").exists() and (path / "documents.json").exists() and (path / "embeddings.npy").exists()


def load_embedding_index(index_dir: str | Path) -> dict | None:
    if not embedding_artifact_exists(index_dir):
        return None
    try:
        import numpy as np
    except ImportError:
        return None

    path = Path(index_dir)
    documents = json.loads((path / "documents.json").read_text(encoding="utf-8"))
    embeddings = np.load(path / "embeddings.npy")
    pmid_to_index = {str(document.get("pmid", "")): index for index, document in enumerate(documents)}
    manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
    return {
        "path": str(path),
        "manifest": manifest,
        "documents": documents,
        "embeddings": embeddings,
        "pmid_to_index": pmid_to_index,
    }


def embedding_coherence(embedding_index: dict | None, pmids: set[str], limit: int = 50) -> float:
    if not embedding_index or len(pmids) < 2:
        return 0.0

    pmid_to_index = embedding_index.get("pmid_to_index", {})
    row_indices = [pmid_to_index[pmid] for pmid in sorted(pmids) if pmid in pmid_to_index][:limit]
    if len(row_indices) < 2:
        return 0.0

    vectors = embedding_index["embeddings"][row_indices]
    similarity = vectors @ vectors.T
    pair_count = len(row_indices) * (len(row_indices) - 1)
    if pair_count <= 0:
        return 0.0
    pair_sum = float(similarity.sum() - similarity.trace())
    return round(max(pair_sum / pair_count, 0.0), 3)


def _result_from_document(
    pmid: str,
    doc: dict,
    score: float,
    query_tokens: list[str],
    reasons: list[str],
    semantic_similarity: float,
) -> dict:
    return {
        "pmid": pmid,
        "score": round(score, 4),
        "title": doc.get("title", ""),
        "publication_year": doc.get("publication_year"),
        "journal": doc.get("journal", ""),
        "publication_types": doc.get("publication_types", []),
        "ranking_reasons": reasons,
        "mesh_terms": doc.get("mesh_terms", [])[:10],
        "authors": doc.get("authors", [])[:5],
        "doi": doc.get("article_ids", {}).get("doi", ""),
        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        "snippet": snippet(doc.get("abstract", ""), query_tokens),
        "semantic_similarity": round(semantic_similarity, 4),
    }


def _document_tokens(document: dict) -> set[str]:
    values = [
        document.get("title", ""),
        " ".join(document.get("mesh_terms", [])),
        " ".join(document.get("publication_types", [])),
    ]
    return set(TOKEN_RE.findall(" ".join(values).lower()))


def query_seed_indices_from_embedding_docs(embedding_index: dict, query_tokens: list[str], seed_limit: int) -> list[int]:
    query_set = {token.lower() for token in query_tokens}
    if not query_set:
        return []
    scored = []
    for row_index, document in enumerate(embedding_index.get("documents", [])):
        overlap = len(query_set & _document_tokens(document))
        if overlap:
            scored.append((overlap, row_index))
    scored.sort(reverse=True)
    return [row_index for _, row_index in scored[:seed_limit]]


def semantic_neighbor_results(
    index: dict,
    embedding_index: dict | None,
    seed_results: list[dict],
    query_tokens: list[str],
    topn: int,
    options: dict,
    seed_limit: int = 25,
    min_similarity: float = 0.6,
) -> list[dict]:
    if not embedding_index or not seed_results or topn <= 0:
        return []

    try:
        import numpy as np
    except ImportError:
        return []

    pmid_to_index = embedding_index.get("pmid_to_index", {})
    seed_pmids = [str(result.get("pmid", "")) for result in seed_results[:seed_limit]]
    seed_indices = [pmid_to_index[pmid] for pmid in seed_pmids if pmid in pmid_to_index]
    if not seed_indices:
        seed_indices = query_seed_indices_from_embedding_docs(embedding_index, query_tokens, seed_limit)
    if not seed_indices:
        return []

    embeddings = embedding_index["embeddings"]
    centroid = np.asarray(embeddings[seed_indices].mean(axis=0), dtype="float32")
    norm = np.linalg.norm(centroid)
    if not norm:
        return []
    centroid = centroid / norm
    scores = embeddings @ centroid
    ranked_indices = np.argsort(-scores)
    documents = index["documents"]
    existing_pmids = {str(result.get("pmid", "")) for result in seed_results}
    year_min = options.get("year_min")
    year_max = options.get("year_max")

    results = []
    for row_index in ranked_indices:
        row_index = int(row_index)
        document_ref = embedding_index["documents"][row_index]
        pmid = str(document_ref.get("pmid", ""))
        if not pmid or pmid in existing_pmids:
            continue
        similarity = float(scores[row_index])
        if similarity < min_similarity:
            break
        doc = documents.get(pmid)
        if not doc:
            continue
        year = doc.get("publication_year")
        if year_min and year and year < year_min:
            continue
        if year_max and year and year > year_max:
            continue
        multiplier, multiplier_reasons = metadata_multiplier(doc, query_tokens, options)
        if multiplier == 0:
            continue
        reasons = ["pubmedbert_semantic_neighbor", *multiplier_reasons]
        results.append(_result_from_document(pmid, doc, similarity * multiplier, query_tokens, reasons, similarity))
        if len(results) >= topn:
            break
    return results
