# Makefile for Protein Interaction Visualizer
# =============================================
# Convenience targets for development and deployment.

.PHONY: help install test clean docker demo notebook lint format

# Default: show help
help:
	@echo "Protein Interaction Visualizer — Make commands"
	@echo ""
	@echo "Setup:"
	@echo "  make install       Create venv and install dependencies"
	@echo "  make dev-install   Install in editable mode with dev deps"
	@echo ""
	@echo "Run:"
	@echo "  make demo          Run synthetic demo pipeline"
	@echo "  make real          Fetch real STRING data and run full pipeline"
	@echo "  make notebook      Start Jupyter Lab"
	@echo ""
	@echo "Testing & Linting:"
	@echo "  make test          Run pytest with coverage"
	@echo "  make lint          Run black, isort, flake8, mypy"
	@echo "  make format        Auto-format code with black+isort"
	@echo ""
	@echo "Docker:"
	@echo "  make docker-build  Build Docker image"
	@echo "  make docker-run    Run container interactively"
	@echo "  make docker-shell  Open bash shell in container"
	@echo ""
	@echo "Cleaning:"
	@echo "  make clean         Remove build artifacts, cache files"
	@echo "  make clean-data    Remove generated data and models"
	@echo "  make clean-all     Deep clean (including venv)"

# --- Environment Setup ---

install:
	@echo "Creating virtual environment..."
	python3 -m venv venv
	@echo "Activate with: source venv/bin/activate"
	@echo "Then run: make dev-install"

dev-install:
	@echo "Installing dependencies..."
	pip install -r requirements.txt
	pip install -e ".[dev]"
	@echo "✅ Done. Run 'make demo' to test."

# --- Execution ---

demo:
	@echo "Running synthetic demo pipeline (no API key needed)..."
	python demo.py

real:
	@echo "Fetching real STRING data and running full pipeline..."
	@echo "Make sure STRING_API_KEY is set in config.yaml or .env"
	python run_pipeline.py

notebook:
	@echo "Starting Jupyter Lab..."
	jupyter lab

# --- Testing & Linting ---

test:
	@echo "Running pytest with coverage..."
	pytest tests/ -v --cov=src --cov-report=term-missing --cov-report=html

lint:
	@echo "Running linters..."
	black --check src/ tests/
	isort --check-only src/ tests/
	flake8 src/ tests/
	mypy src/ 2>/dev/null || echo "Mypy not configured; skipping"

format:
	@echo "Formatting code..."
	black src/ tests/
	isort src/ tests/

# --- Docker ---

docker-build:
	@echo "Building Docker image..."
	docker build -t ppi-visualizer .

docker-run:
	@echo "Running container (mounts ./data)..."
	docker run -p 8888:8888 -v $$(pwd)/data:/app/data ppi-visualizer

docker-shell:
	@echo "Opening shell in container..."
	docker run -it --entrypoint /bin/bash ppi-visualizer

# --- Cleaning ---

clean:
	@echo "Cleaning Python artifacts..."
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	rm -rf .mypy_cache .pytest_cache .coverage htmlcov
	@echo "✅ Clean"

clean-data:
	@echo "WARNING: Removing all generated data and models"
	rm -rf data/processed/*.csv data/processed/*.json models/*.joblib results/figures/*.html results/predictions/*.csv logs/*.log
	@echo "✅ Data cleaned (directories kept)"

clean-all: clean clean-data
	@echo "Removing virtual environment..."
	rm -rf venv
	@echo "✅ Full clean complete"

# --- Utilities ---

# Show current config
config:
	@cat config.yaml

# Show what data sources are available
sources:
	@echo "Available data sources:"
	@echo "  - STRING (https://string-db.org)"
	@echo "  - BioGRID (https://thebiogrid.org)"
	@echo "Get STRING API key: https://string-db.org/cgi/access?section=api"

# Install pre-commit hooks
pre-commit:
	pre-commit install
	@echo "Pre-commit hooks installed. They will run on git commit."

# Generate requirements.lock from requirements.in (if using pip-tools)
lock:
	pip-compile requirements.in -o requirements.txt 2>/dev/null || echo "No requirements.in or pip-tools not installed. Skipping."
