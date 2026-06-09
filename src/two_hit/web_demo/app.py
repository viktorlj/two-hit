"""FastAPI application factory for the lightweight two-hit demo.

No idle watchdog, no PDF, no JSON API, no result store — kept small enough to
run on a modest Railway container.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .. import __version__

WEB_DIR = Path(__file__).parent
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"


def create_app() -> FastAPI:
    app = FastAPI(
        title="two-hit (demo)",
        description="Lightweight demo of integrated mutation + copy number analysis",
        version=__version__,
    )
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    templates = Jinja2Templates(directory=str(TEMPLATE_DIR))
    templates.env.globals["version"] = __version__
    app.state.templates = templates

    from .routes import router

    app.include_router(router)
    return app


app = create_app()
