"""Tests for src.network_analysis.graph_analyzer.

Tests:
    - GraphAnalyzer initialization
    - Metric computation (density, degree, clustering, etc.)
    - Community detection
    - Hub identification

Run with: pytest tests/test_network_analysis.py -v
"""

import pytest
import networkx as nx
import numpy as np

from src.network_analysis.graph_analyzer import GraphAnalyzer, load_network_from_csv
import pandas as pd
import tempfile
from pathlib import Path


class TestGraphAnalyzerInit:
    """GraphAnalyzer construction and basic properties."""

    def test_init_with_networkx_graph(self):
        """Initialize from NetworkX graph."""
        G = nx.karate_club_graph()
        analyzer = GraphAnalyzer(G)
        assert analyzer.G.number_of_nodes() == 34
        assert analyzer.G.number_of_edges() == 78

    def test_init_with_directed_graph(self):
        """Directed graphs are converted to undirected by default."""
        G = nx.DiGraph()
        G.add_edges_from([(1, 2), (2, 3)])
        analyzer = GraphAnalyzer(G, directed=False)
        assert not analyzer.G.is_directed()

    def test_empty_graph_warning(self):
        """Empty graph produces sensible metrics."""
        G = nx.Graph()
        analyzer = GraphAnalyzer(G)
        metrics = analyzer.compute_metrics()
        assert metrics["num_nodes"] == 0
        assert metrics["num_edges"] == 0


class TestNetworkMetrics:
    """Network metric computation."""

    @pytest.fixture
    def small_graph(self):
        """Create a small test graph (Karate Club)."""
        G = nx.karate_club_graph()
        return GraphAnalyzer(G)

    def test_density(self, small_graph):
        """Density is between 0 and 1."""
        metrics = small_graph.compute_metrics()
        assert 0 <= metrics["density"] <= 1

    def test_degree_metrics(self, small_graph):
        """Degree statistics are computed."""
        metrics = small_graph.compute_metrics()
        assert metrics["degree_mean"] > 0
        assert metrics["degree_max"] >= metrics["degree_mean"]
        assert metrics["degree_median"] > 0

    def test_clustering(self, small_graph):
        """Clustering coefficients are computed."""
        metrics = small_graph.compute_metrics()
        assert 0 <= metrics["transitivity"] <= 1
        assert 0 <= metrics["avg_clustering_coefficient"] <= 1

    def test_connected_components(self, small_graph):
        """Connected components are counted."""
        metrics = small_graph.compute_metrics()
        assert metrics["num_connected_components"] >= 1
        assert metrics["largest_component_size"] > 0

    def test_path_metrics(self, small_graph):
        """Diameter and shortest path are positive."""
        metrics = small_graph.compute_metrics()
        if metrics["diameter"] is not None:
            assert metrics["diameter"] >= 0
        if metrics["avg_shortest_path_length"] is not None:
            assert metrics["avg_shortest_path_length"] >= 0

    def test_centrality_computed(self, small_graph):
        """Centrality measures present if enabled in config."""
        metrics = small_graph.compute_metrics()
        # Degree centrality always on
        assert "degree_centrality" in metrics
        # Check we got actual centrality values
        assert len(metrics["degree_centrality"]) == small_graph.G.number_of_nodes()


class TestCommunityDetection:
    """Community detection algorithms."""

    @pytest.fixture
    def modular_graph(self):
        """Create a graph with known community structure (2 communities)."""
        G = nx.Graph()
        # Community 1: clique of 5
        G.add_edges_from([(i, j) for i in range(5) for j in range(i + 1, 5)])
        # Community 2: clique of 5
        G.add_edges_from([(i + 5, j + 5) for i in range(5) for j in range(i + 1, 5)])
        # Weak bridge
        G.add_edge(0, 5)
        return GraphAnalyzer(G)

    def test_louvain_detection(self, modular_graph):
        """Louvain detects communities."""
        communities = modular_graph.detect_communities(algorithm="louvain")
        assert len(communities) >= 2  # At least 2 groups
        # Check all nodes assigned
        all_nodes = set(modular_graph.G.nodes())
        assigned = set()
        for nodes in communities.values():
            assigned.update(nodes)
        assert assigned == all_nodes

    def test_label_propagation(self, modular_graph):
        """Label propagation runs without error."""
        communities = modular_graph.detect_communities(algorithm="label_propagation")
        assert isinstance(communities, dict)

    def test_girvan_newman_slow(self, modular_graph):
        """Girvan-Newman works (but returns top-level split for small graphs)."""
        communities = modular_graph.detect_communities(algorithm="girvan_newman")
        assert isinstance(communities, dict)


class TestHubIdentification:
    """Hub protein identification."""

    def test_find_hubs_returns_dataframe(self):
        """find_hubs returns sorted DataFrame."""
        G = nx.barabasi_albert_graph(100, 3, seed=42)
        analyzer = GraphAnalyzer(G)
        hubs = analyzer.find_hubs(top_n=10)
        assert isinstance(hubs, pd.DataFrame)
        assert len(hubs) == 10
        # Sorted descending by centrality
        assert hubs["centrality_score"].iloc[0] >= hubs["centrality_score"].iloc[-1]
        # Has required columns
        assert {"protein_id", "centrality_score", "degree", "community_id"}.issubset(hubs.columns)


class TestNetworkIO:
    """Loading networks from files."""

    def test_load_from_csv(self, tmp_path):
        """CSV loader creates graph correctly."""
        csv_path = tmp_path / "test.csv"
        csv_path.write_text("protein1,protein2,combined_score\nA,B,900\nB,C,800\n")
        G = load_network_from_csv(csv_path)
        assert G.number_of_nodes() == 3
        assert G.number_of_edges() == 2
        data = G["A"]["B"]
        assert data["combined_score"] == 900


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
