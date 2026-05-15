"""Pytest configuration for the test suite.

Ensures src/ is on the Python path so imports like 'from src.X' work
without needing pip install -e . first.
"""

import sys
from pathlib import Path

# Add the src/ directory to the path so 'from src.X' imports work
# This is needed for CI where we don't do 'pip install -e .'
PROJECT_ROOT = Path(__file__).parent.parent
SRC_DIR = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))
