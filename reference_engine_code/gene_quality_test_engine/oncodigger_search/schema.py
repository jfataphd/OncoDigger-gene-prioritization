from dataclasses import dataclass, field


@dataclass
class SearchResult:
    pmid: str
    score: float
    title: str
    publication_year: int | None = None
    journal: str = ""
    publication_types: list[str] = field(default_factory=list)
    mesh_terms: list[str] = field(default_factory=list)
    authors: list[dict] = field(default_factory=list)
    doi: str = ""
    pubmed_url: str = ""
    snippet: str = ""
    ranking_reasons: list[str] = field(default_factory=list)


def result_to_dict(result: SearchResult) -> dict:
    return {
        "pmid": result.pmid,
        "score": result.score,
        "title": result.title,
        "publication_year": result.publication_year,
        "journal": result.journal,
        "publication_types": result.publication_types,
        "mesh_terms": result.mesh_terms,
        "authors": result.authors,
        "doi": result.doi,
        "pubmed_url": result.pubmed_url,
        "snippet": result.snippet,
        "ranking_reasons": result.ranking_reasons,
    }
