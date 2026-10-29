from collections import Counter
import re


ARTICLE_TYPE_OPTIONS = [
    "Research articles",
    "Clinical trials",
    "Case reports",
    "Reviews and evidence syntheses",
]

DEFAULT_ARTICLE_TYPES = list(ARTICLE_TYPE_OPTIONS)

REVIEW_TITLE_RE = re.compile(
    r"\b(review|systematic review|meta-analysis|meta analysis|literature review|overview|state of the art)\b",
    re.IGNORECASE,
)


def classify_article_types(doc: dict) -> set[str]:
    cached = doc.get("_oncodigger_article_types")
    if cached is not None:
        return set(cached)
    publication_types = {publication_type.lower() for publication_type in doc.get("publication_types", [])}
    title = doc.get("title", "")
    categories = set()

    if publication_types & {"clinical trial", "controlled clinical trial", "randomized controlled trial"}:
        categories.add("Clinical trials")
    if any("clinical trial" in publication_type for publication_type in publication_types):
        categories.add("Clinical trials")
    if "case reports" in publication_types:
        categories.add("Case reports")
    if publication_types & {"systematic review", "meta-analysis"}:
        categories.add("Reviews and evidence syntheses")
    elif "review" in publication_types and REVIEW_TITLE_RE.search(title):
        categories.add("Reviews and evidence syntheses")

    research_signals = {
        "journal article",
        "comparative study",
        "evaluation study",
        "validation study",
        "observational study",
        "multicenter study",
    }
    if publication_types & research_signals:
        categories.add("Research articles")
    if not categories:
        categories.add("Research articles")
    doc["_oncodigger_article_types"] = tuple(sorted(categories))
    return categories


def filter_results_by_article_type(index, results: list[dict], selected_article_types: list[str]) -> list[dict]:
    if not selected_article_types:
        return []
    if set(selected_article_types) == set(ARTICLE_TYPE_OPTIONS):
        return results
    selected = set(selected_article_types)
    return [result for result in results if classify_article_types(result_document(index, result)) & selected]


def result_document(index, result: dict) -> dict:
    if "publication_types" in result and "title" in result:
        return result
    # Offline callers often supply only PMIDs, unlike the Streamlit BM25 results.
    from .discovery import _get_doc
    return {**_get_doc(index, result["pmid"]), **result}


def count_article_types(index, results: list[dict]) -> Counter:
    counts = Counter()
    for result in results:
        for article_type in classify_article_types(result_document(index, result)):
            if article_type in ARTICLE_TYPE_OPTIONS:
                counts[article_type] += 1
    return counts


def result_manuscript_types(index, result: dict) -> str:
    article_types = [article_type for article_type in ARTICLE_TYPE_OPTIONS if article_type in classify_article_types(result_document(index, result))]
    return ", ".join(article_types)
