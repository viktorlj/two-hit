"""FastAPI application for two-hit web interface."""

from __future__ import annotations

import asyncio
import logging
import os
import time
import uuid
from collections import OrderedDict
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .. import __version__

logger = logging.getLogger(__name__)

WEB_DIR = Path(__file__).parent
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"
MAX_RESULTS = 500

# Idle timeout: exit with non-zero code after this many seconds without a
# request, so the container platform (Railway, Cloud Run) restarts the process
# with a clean memory slate.  Set to 0 to disable.
IDLE_TIMEOUT = int(os.environ.get("IDLE_TIMEOUT", "1800"))  # 30 min default


class ResultStore:
    """Simple in-memory result store with FIFO eviction."""

    def __init__(self, maxsize: int = MAX_RESULTS):
        self._store: OrderedDict[str, dict] = OrderedDict()
        self._maxsize = maxsize

    def put(self, data: dict) -> str:
        rid = uuid.uuid4().hex[:12]
        if len(self._store) >= self._maxsize:
            self._store.popitem(last=False)
        self._store[rid] = data
        return rid

    def get(self, rid: str) -> dict | None:
        return self._store.get(rid)


_last_request_time: float = time.monotonic()


async def _idle_watchdog() -> None:
    """Background task that exits the process after IDLE_TIMEOUT seconds of inactivity."""
    while True:
        await asyncio.sleep(60)  # check every minute
        idle = time.monotonic() - _last_request_time
        if idle > IDLE_TIMEOUT:
            logger.info("Idle for %.0fs (limit %ds), exiting for restart", idle, IDLE_TIMEOUT)
            os._exit(1)  # non-zero so Railway restarts the container


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Start the idle watchdog on app startup."""
    task = None
    if IDLE_TIMEOUT > 0:
        task = asyncio.create_task(_idle_watchdog())
        logger.info("Idle watchdog started (timeout=%ds)", IDLE_TIMEOUT)
    yield
    if task:
        task.cancel()


def create_app() -> FastAPI:
    """Create and configure the FastAPI application."""
    app = FastAPI(
        title="two-hit",
        description="Integrated mutation + copy number biallelic inactivation analysis",
        version=__version__,
        lifespan=lifespan,
    )

    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

    templates = Jinja2Templates(directory=str(TEMPLATE_DIR))
    templates.env.globals["version"] = __version__

    app.state.templates = templates
    app.state.results = ResultStore()

    @app.middleware("http")
    async def track_activity(request: Request, call_next):
        global _last_request_time
        _last_request_time = time.monotonic()
        return await call_next(request)

    from .routes.analyze import router as analyze_router

    app.include_router(analyze_router)

    return app


app = create_app()
