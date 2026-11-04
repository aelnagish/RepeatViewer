"""Find tips in a tree by name or trait value (no Qt)."""

from __future__ import annotations

import re


def split_terms(query: str) -> list[str]:
    """Words separated by spaces, commas or semicolons; a tip matches any of them."""
    return [t for t in re.split(r"[\s,;]+", query or "") if t]


def find_tips(
    tips: list[str],
    query: str,
    *,
    traits: dict[str, dict[str, str]] | None = None,
    in_names: bool = True,
    in_traits: bool = True,
    case_sensitive: bool = False,
) -> list[str]:
    """Tips (in the order given) whose name or trait values contain a search word."""
    terms = split_terms(query)
    if not terms:
        return []
    if not case_sensitive:
        terms = [t.lower() for t in terms]

    def hay(text: str) -> str:
        return text if case_sensitive else text.lower()

    hits: list[str] = []
    for tip in tips:
        texts: list[str] = []
        if in_names:
            texts.append(tip)
        if in_traits and traits:
            texts.extend(str(v) for v in traits.get(tip, {}).values())
        if any(term in hay(text) for term in terms for text in texts):
            hits.append(tip)
    return hits


def step_index(current: int, count: int, direction: int) -> int:
    """Next (+1) or previous (-1) hit index, wrapping; ``current`` -1 = nothing shown yet."""
    if count <= 0:
        return -1
    if current < 0:
        return 0 if direction >= 0 else count - 1
    return (current + (1 if direction >= 0 else -1)) % count
