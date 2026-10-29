# proximity.py

Standalone proximity module written by the author in the application checkout (`src/oncodigger/search/proximity.py`,
uncommitted at the time of copying, 2026-10-14) and copied here unchanged. It computes the per-paper proximity 1/(1+d), the mean over supporting papers
and the adjusted relevance (relevance x (1 + w x P), w = 3 by default). It is used by `scripts/analysis/proximity_boost_run.py` for the post hoc
analysis in Methods 2.8 / Results 5.11 and is engine-independent.

sha256 of the copied file: `d8103d24c3cf3034898b6b0dbdd104995b026829263ba031e9285c33206332bc`
