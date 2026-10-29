from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Mapping


@dataclass(frozen=True)
class ProximityResult:
    proximity: float
    adjusted_relevance: float


def term_spans(tokens: list[str], terms: Iterable[str]) -> list[tuple[int, int]]:
    spans = []
    for term in terms:
        term_tokens = term.split()
        if not term_tokens:
            continue
        width = len(term_tokens)
        for start in range(len(tokens) - width + 1):
            if tokens[start:start + width] == term_tokens:
                spans.append((start, start + width))
    return spans


def span_distance(left: tuple[int, int], right: tuple[int, int]) -> int:
    if left[1] <= right[0]:
        return right[0] - left[1] + 1
    if right[1] <= left[0]:
        return left[0] - right[1] + 1
    return 0


def document_proximity(
    paper: Mapping[str, object],
    query_terms: set[str],
    entity_terms: set[str],
    normalize_segment: Callable[[str], str],
) -> float:
    """Best title/abstract proximity, with the two fields treated separately."""
    best = None
    for field in ("title", "abstract"):
        tokens = normalize_segment(str(paper.get(field, "") or "")).split()
        query_spans = term_spans(tokens, query_terms)
        entity_spans = term_spans(tokens, entity_terms)
        if not query_spans or not entity_spans:
            continue
        distance = min(
            span_distance(query_span, entity_span)
            for query_span in query_spans
            for entity_span in entity_spans
        )
        best = distance if best is None else min(best, distance)
    return 0.0 if best is None else 1.0 / (1.0 + best)


def score_supporting_papers(
    supporting_papers: Iterable[Mapping[str, object]],
    query_terms: set[str],
    entity_terms_by_pmid: Mapping[str, set[str]],
    normalize_segment: Callable[[str], str],
    original_relevance: float,
    weight: float = 3.0,
) -> ProximityResult:
    """Apply the frozen continuous definition without changing evidence selection."""
    values = []
    for paper in supporting_papers:
        pmid = str(paper.get("pmid", ""))
        values.append(
            document_proximity(
                paper,
                query_terms,
                entity_terms_by_pmid.get(pmid, set()),
                normalize_segment,
            )
        )
    proximity = sum(values) / len(values) if values else 0.0
    return ProximityResult(
        proximity=proximity,
        adjusted_relevance=original_relevance * (1.0 + weight * proximity),
    )
