FROM python:3.12-slim

WORKDIR /app

# Install uv for fast dependency resolution
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Copy project files
COPY pyproject.toml README.md ./
COPY src/ src/
COPY demo/ demo/
COPY data/ data/

# Install the package (system-wide, no venv needed in container)
RUN uv pip install --system --no-cache .

# Railway sets PORT env var; default to 8080
ENV PORT=8080
EXPOSE 8080

CMD exec uvicorn two_hit.web.app:app --host 0.0.0.0 --port $PORT
