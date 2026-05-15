#!/usr/bin/env python3
"""
setup.py — Installable package configuration for Protein Interaction Visualizer

This makes the project importable as a Python package:
    pip install -e .

Then you can:
    from src.data_acquisition.string_fetcher import StringFetcher

For more details see pyproject.toml (PEP 621 metadata).
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="protein-interaction-visualizer",
    version="0.1.0",
    author="Your Name",  # TODO: Replace with your name
    author_email="your.email@example.com",  # TODO: Replace
    description="A toolkit for fetching, visualizing, and predicting protein-protein interactions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/protein-interaction-visualizer",  # TODO: Replace
    project_urls={
        "Bug Tracker": "https://github.com/yourusername/protein-interaction-visualizer/issues",
        "Documentation": "https://github.com/yourusername/protein-interaction-visualizer/docs",
        "Source Code": "https://github.com/yourusername/protein-interaction-visualizer",
    },
    classifiers=[
        "Development Status :: 3 - Alpha",  # Change to "4 - Beta" when stable
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.9, <3.13",
    install_requires=[
        "numpy>=1.24.0",
        "pandas>=2.0.0",
        "scipy>=1.10.0",
        "networkx>=3.0",
        "scikit-learn>=1.3.0",
        "plotly>=5.14.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "requests>=2.30.0",
        "python-dotenv>=1.0.0",
        "pyyaml>=6.0",
        "joblib>=1.3.0",
        "tqdm>=4.65.0",
        "biopython>=1.81",  # Optional but useful for sequence data
        "colorlog>=6.7.0",
        "python-dateutil>=2.8.0",
    ],
    extras_require={
        # Optional dependencies for extra features
        "dev": [  # Development tools (testing, linting, formatting)
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.4.0",
        ],
        "notebooks": [  # Jupyter notebooks support
            "jupyterlab>=4.0.0",
            "ipykernel>=6.25.0",
            "nbconvert>=7.7.0",
        ],
        "graphviz": [  # Fast graph layout rendering
            "pygraphviz; platform_system!='Windows'",  # Requires system Graphviz
        ],
        "deep-learning": [  # Graph Neural Networks
            "torch>=2.0.0",
            "torch-geometric>=2.3.0",
        ],
        "bio": [  # Advanced bioinformatics
            "goatools>=0.1.9",  # Gene Ontology analysis
            "gseapy>=1.0.6",  # Gene set enrichment
        ],
    },
    entry_points={
        # Command-line scripts installed to user's PATH
        "console_scripts": [
            "ppi-fetch=src.data_acquisition.string_fetcher:main",
            "ppi-viz=src.visualization.network_viz:main",
            "ppi-predict=src.prediction.interaction_predictor:main",
            "ppi-analyze=src.network_analysis.graph_analyzer:main",
        ],
    },
    include_package_data=True,
    package_data={
        # Include config files, example data
        "protein_interaction_visualizer": [
            "config.yaml",
            ".env.example",
        ],
    },
)
