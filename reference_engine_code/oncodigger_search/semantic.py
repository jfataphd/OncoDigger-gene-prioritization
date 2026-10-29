from collections import Counter
import math


def semantic_neighborhood_score(
    bridge_terms: str,
    support_percent: float,
    corpus_percent: float,
    supporting_papers: int,
    embedding_coherence: float = 0.0,
) -> float:
    """Estimate concept-neighborhood strength, using precomputed embeddings when available."""
    bridge_count = len([term for term in bridge_terms.split(", ") if term])
    enrichment_gap = max(support_percent - corpus_percent, 0)
    support_weight = math.log1p(max(supporting_papers, 0))
    embedding_weight = max(embedding_coherence, 0.0) * 3.0
    return round((bridge_count * 0.4) + (enrichment_gap * 0.08) + support_weight + embedding_weight, 3)


def semantic_signal(semantic_score: float, relationship_signal: str) -> str:
    if semantic_score >= 8:
        return "Strong concept neighborhood"
    if relationship_signal == "Enriched relationship" and semantic_score >= 5:
        return "Enriched concept neighborhood"
    if semantic_score >= 3:
        return "Moderate concept neighborhood"
    return "Sparse concept neighborhood"


def shared_bridge_terms(rows: list[dict], limit: int = 20) -> list[dict]:
    counts = Counter()
    examples = {}
    for row in rows:
        label = row.get("Term") or row.get("Gene") or row.get("Drug") or row.get("Phytochemical") or row.get("Compound") or ""
        for term in str(row.get("Bridge terms", "")).split(", "):
            if not term:
                continue
            counts[term] += 1
            examples.setdefault(term, []).append(label)
    return [
        {
            "Bridge term": term,
            "Relationship count": count,
            "Example relationships": ", ".join(examples.get(term, [])[:5]),
        }
        for term, count in counts.most_common(limit)
    ]
