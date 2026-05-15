# Contributing to Protein Interaction Visualizer

First off, thank you for considering contributing! 🎉

This document provides guidelines for developers who want to improve the project.

---

## Code of Conduct

Be respectful, inclusive, and constructive. Harassment or toxic behavior will not be tolerated.

---

## How to Contribute

### Reporting Bugs
- **Search existing issues** first to avoid duplicates
- Use the bug report template
- Include:
  - OS, Python version, package versions (`pip freeze`)
  - Steps to reproduce (minimal example)
  - Expected vs actual behavior
  - Error messages and stack traces
  - Log files (`logs/`)

### Suggesting Features
- Open a "Feature request" issue
- Describe the use case and expected API
- Discuss design before large PRs

### Pull Requests

1. **Fork** the repository and create a feature branch:
   ```bash
   git checkout -b feature/my-cool-additions
   ```

2. **Set up development environment:**
   ```bash
   pip install -e ".[dev]"  # includes pytest, black, isort, mypy
   pre-commit install        # optional: run checks on commit
   ```

3. **Write clean, commented code:**
   - Follow PEP 8 (line length 100)
   - Add docstrings to all public functions/classes (Google or NumPy style)
   - Type hints where practical

4. **Add tests:** Aim for at least 80% coverage on new code.
   ```bash
   pytest tests/test_my_module.py -v
   ```

5. **Run linters:**
   ```bash
   black src/ tests/
   isort src/ tests/
   flake8 src/ tests/
   mypy src/  # optional but nice
   ```

6. **Commit with clear messages:**
   ```
   feat: add BioGRID API v3.5 support
   fix: handle missing community labels in visualization
   docs: update installation instructions
   ```

7. **Push and open PR** against `main` branch.
   - Include before/after screenshots if UI change
   - Reference related issues (e.g., "Fixes #123")

8. **CI will run automatically.** Address any failing checks.

---

## Development Workflow

### Directory Layout
- New modules go in `src/` (as subpackages)
- Tests mirror structure: `tests/test_<module>.py`
- Docs: docs/ (Markdown or Sphinx)
- Notebooks: notebooks/ (for demos, tutorials)

### Adding a New Data Fetcher

1. Create `src/data_acquisition/myfetcher.py`
2. Subclass `BaseFetcher` if appropriate, or follow same interface:
   - `__init__(self, param1=default)`
   - `fetch_interactions(self) -> pd.DataFrame`
   - Standardized columns: `["protein1", "protein2", "combined_score", "source"]`
3. Add to `src/data_acquisition/__init__.py`
4. Write unit tests in `tests/test_data_acquisition.py`
5. Update `config.yaml` schema and `README.md`

### Adding a New Visualization Type

1. Add method to `NetworkVisualizer` in `src/visualization/network_viz.py`
2. Document parameters in docstring
3. Add test in `tests/test_visualization.py`
4. Optionally add to notebook demo

### Adding a New ML Model

1. Extend `InteractionPredictor.train()` to accept new model type
2. Implement scikit-learn compatible API: `fit`, `predict`, `predict_proba`
3. Add hyperparameter config section in `config.yaml`
4. Add model-specific evaluation (e.g., learning curves)
5. Document when to use it vs default

---

## Release Process

Versioning follows SemVer: `MAJOR.MINOR.PATCH`

**When to bump:**
- MAJOR: breaking API change
- MINOR: new feature, backward-compatible
- PATCH: bug fix

Steps:
1. Update `CHANGELOG.md` with changes
2. Update version in `pyproject.toml` and `src/__init__.py`
3. Create GitHub release with release notes
4. Publish to PyPI (if applicable):
   ```bash
   python -m build
   twine upload dist/*
   ```

---

## Style Guide

### Python Style
- 4 spaces indent
- Max line length: 100 characters (black default)
- f-strings for formatting
- Use pathlib, not os.path
- Prefer data classes over tuples for multi-field data

### Docstring Format (Google style example)
```python
def plot_interactive_graph(
    title: str,
    output_path: Path | None = None,
) -> go.Figure:
    """Create an interactive Plotly network figure.

    Args:
        title: Plot title displayed at top.
        output_path: If given, HTML saved here. Directory created automatically.
    Returns:
        Plotly Figure object (also saved if output_path provided).
    """
```

### Imports (isort order)
```python
# Standard library
import json
import logging
from pathlib import Path
from typing import List

# Third-party
import networkx as nx
import pandas as pd
import plotly.express as px

# Local (project)
from src.utils.config import load_config
```

### Logging
Use module-level logger:
```python
from src.utils.logging import get_logger
logger = get_logger(__name__)

logger.info("Network loaded: %d nodes", G.number_of_nodes())
logger.debug("Detailed %s", expensive_var)
logger.warning("Score threshold high: %d", threshold)
logger.error("Failed to fetch: %s", error)
```

---

## Dependency Management

- **Core dependencies** pinned in `requirements.txt`.  
  Update after testing: `pip-compile requirements.in` if using `pip-tools`.
- **Optional dependencies** noted in `extras_require` in `pyproject.toml`.
- **Development dependencies** (`dev` extra): pytest, black, isort, mypy, pre-commit.

Add new packages by:
1. Adding to appropriate section in `pyproject.toml`
2. Running `pip install -e ".[dev]"` locally
3. Updating `requirements.txt` if you maintain that lockfile

---

## Performance Tips

- For analyses > 10k nodes, set `compute_betweenness=False` in config
- Use `layout_algorithm: "kk"` for faster layout on moderate graphs
- For ML training, subsample edges if graph too large (e.g., `max_edges=50000`)
- Cache intermediate results (CSV) rather than re-fetching from API

---

## Questions?

Open an issue or reach out via email.

Happy coding! 🎯
