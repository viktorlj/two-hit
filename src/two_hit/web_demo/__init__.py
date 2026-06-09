"""Lightweight two-hit demo web app (Railway-friendly).

Drops plotly and fpdf2; restricts analysis to a curated gene panel and
bounded uploads. The full app lives in two_hit.web.
"""

from __future__ import annotations


def create_app():
    """Lazily build the demo FastAPI app (avoids import cost at package import)."""
    from .app import create_app as _create

    return _create()
