import math
import sqlite3


def enrichment_score(
    weighted_score: float,
    supporting_papers: int,
    relevant_papers: int,
    background_papers: int,
    corpus_papers: int,
) -> float:
    """Rank query-specific associations against the full corpus background."""
    if supporting_papers <= 0 or relevant_papers <= 0:
        return 0.0
    corpus_papers = max(corpus_papers, 1)
    background_papers = max(background_papers, supporting_papers, 1)
    relevant_fraction = supporting_papers / relevant_papers
    background_fraction = background_papers / corpus_papers
    enrichment = relevant_fraction / max(background_fraction, 1 / corpus_papers)
    support_boost = math.log1p(supporting_papers)
    return weighted_score * math.log1p(enrichment) * support_boost


def single_term_background_papers(index, term: str) -> int:
    parts = term.split()
    if len(parts) != 1:
        return 0
    if isinstance(index, sqlite3.Connection):
        row = index.execute("SELECT COUNT(*) FROM postings WHERE token = ?", (parts[0],)).fetchone()
        return row[0] if row else 0
    return len(index.get("postings", {}).get(parts[0], []))


def entity_background_papers(index: dict, terms: list[str]) -> int:
    background_counts = [single_term_background_papers(index, term) for term in terms]
    background_counts = [count for count in background_counts if count > 0]
    return min(background_counts) if background_counts else 0
