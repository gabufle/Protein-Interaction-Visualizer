"""Graph analysis utilities for protein interaction networks.

Compute network topology metrics, identify key proteins (hubs),
detect communities, and quantify network properties.

Based on NetworkX graph algorithms and standard network science.

Author: [Your Name]
Date: 2025-05-14
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from community import community_louvain  # python-louvain package
from tqdm import tqdm

from src.utils.config import load_config, get_random_seed
from src.utils.logging import get_logger

logger = get_logger(__name__)


class GraphAnalyzer:
    """
    Comprehensive analysis of a protein-protein interaction network.

    Computes standard network metrics:
        - Basic: nodes, edges, density, connected components
        - Centrality: degree, betweenness, closeness, eigenvector, Katz
        - Clustering: clustering coefficient, transitivity
        - Communities: modularity-based detection (Louvain, Girvan-Newman)
        - Hubs: identifies connector proteins by centrality

    All results are returned as dict for easy serialization.

    Usage:
        import networkx as nx
        from src.network_analysis.graph_analyzer import GraphAnalyzer

        G = nx.from_pandas_edgelist(df, "protein1", "protein2")
        analyzer = GraphAnalyzer(G, directed=False)
        metrics = analyzer.compute_metrics()
        hubs = analyzer.find_hubs(top_n=20)
        analyzer.save_results("results/network_analysis.json")
    """

    def __init__(
        self,
        graph: nx.Graph,
        directed: bool = False,
        config: Optional[Dict] = None,
    ):
        """
        Initialize analyzer.

        Args:
            graph: NetworkX graph (undirected or directed)
            directed: If True, treat as directed graph (else treat as undirected)
            config: Optional config dict overriding defaults
        """
        self.G = graph
        self.directed = directed and graph.is_directed()

        if not self.directed:
            # Ensure undirected for consistency
            self.G = self.G.to_undirected()

        self.config = config or load_config()
        self.random_seed = get_random_seed()

        logger.info(
            f"GraphAnalyzer initialized: "
            f"{self.G.number_of_nodes():,} nodes, {self.G.number_of_edges():,} edges"
        )

    def compute_metrics(self) -> Dict:
        """
        Compute core network metrics.

        Returns:
            Dict of metric_name -> value (int/float)
        """
        logger.info("Computing network metrics...")
        metrics = {}

        # Basic statistics
        metrics["num_nodes"] = self.G.number_of_nodes()
        metrics["num_edges"] = self.G.number_of_edges()
        metrics["density"] = nx.density(self.G)

        # Connected components
        if self.G.is_directed():
            metrics["num_strongly_connected_components"] = nx.number_strongly_connected_components(self.G)
            metrics["num_weakly_connected_components"] = nx.number_weakly_connected_components(self.G)
            largest_cc = max(nx.weakly_connected_components(self.G), key=len)
        else:
            metrics["num_connected_components"] = nx.number_connected_components(self.G)
            largest_cc = max(nx.connected_components(self.G), key=len)

        metrics["largest_component_size"] = len(largest_cc)
        metrics["largest_component_ratio"] = len(largest_cc) / self.G.number_of_nodes()

        # Path lengths (only on largest CC to avoid infinity)
        G_largest = self.G.subgraph(largest_cc).copy()
        if G_largest.number_of_nodes() > 1:
            try:
                metrics["diameter"] = nx.diameter(G_largest)
                metrics["radius"] = nx.radius(G_largest)
                metrics["avg_shortest_path_length"] = nx.average_shortest_path_length(G_largest)
            except Exception as e:
                logger.warning(f"Could not compute some path metrics: {e}")
                metrics["diameter"] = None
                metrics["radius"] = None
                metrics["avg_shortest_path_length"] = None

        # Clustering
        metrics["transitivity"] = nx.transitivity(self.G)  # global clustering coeff
        metrics["avg_clustering_coefficient"] = nx.average_clustering(self.G)

        # Assortativity (degree correlation)
        try:
            metrics["degree_assortativity"] = nx.degree_assortativity_coefficient(self.G)
        except Exception:
            metrics["degree_assortativity"] = None

        # Centrality (expensive; compute selectively)
        centrality_cfg = self.config.get("network", {}).get("centrality", {})

        if centrality_cfg.get("compute_degree", True):
            logger.debug("Computing degree centrality...")
            deg_cent = nx.degree_centrality(self.G)
            metrics["degree_centrality_mean"] = np.mean(list(deg_cent.values()))
            metrics["degree_centrality_max"] = max(deg_cent.values())
            metrics["degree_centrality"] = deg_cent  # full dict

        if centrality_cfg.get("compute_betweenness", False):
            logger.debug("Computing betweenness centrality (this may take a while)...")
            k = centrality_cfg.get("betweenness_k", None)  # approximate; None = exact
            betw_cent = nx.betweenness_centrality(self.G, k=k, seed=self.random_seed)
            metrics["betweenness_centrality_mean"] = np.mean(list(betw_cent.values()))
            metrics["betweenness_centrality_max"] = max(betw_cent.values())
            metrics["betweenness_centrality"] = betw_cent

        if centrality_cfg.get("compute_closeness", True):
            logger.debug("Computing closeness centrality...")
            close_cent = nx.closeness_centrality(self.G)
            metrics["closeness_centrality_mean"] = np.mean(list(close_cent.values()))
            metrics["closeness_centrality_max"] = max(close_cent.values())
            metrics["closeness_centrality"] = close_cent

        if centrality_cfg.get("compute_eigenvector", False):
            logger.debug("Computing eigenvector centrality...")
            try:
                eig_cent = nx.eigenvector_centrality_numpy(self.G, weight="weight")
                metrics["eigenvector_centrality_mean"] = np.mean(list(eig_cent.values()))
                metrics["eigenvector_centrality_max"] = max(eig_cent.values())
                metrics["eigenvector_centrality"] = eig_cent
            except Exception as e:
                logger.warning(f"Eigenvector centrality failed (non-convergent?): {e}")
                metrics["eigenvector_centrality"] = {}

        if centrality_cfg.get("compute_katz", False):
            logger.debug("Computing Katz centrality...")
            try:
                katz_cent = nx.katz_centrality_numpy(self.G, weight="weight")
                metrics["katz_centrality_mean"] = np.mean(list(katz_cent.values()))
                metrics["katz_centrality_max"] = max(katz_cent.values())
                metrics["katz_centrality"] = katz_cent
            except Exception as e:
                logger.warning(f"Katz centrality failed: {e}")
                metrics["katz_centrality"] = {}

        # Degree distribution
        degrees = [d for _, d in self.G.degree()]
        metrics["degree_mean"] = np.mean(degrees)
        metrics["degree_median"] = np.median(degrees)
        metrics["degree_max"] = max(degrees)
        metrics["degree_std"] = np.std(degrees)

        # Power law fit (scale-free property) — estimate exponent
        try:
            from scipy import stats

            # Fit power law to degree distribution (excluding degree=0)
            nonzero_degrees = [d for d in degrees if d > 0]
            if len(nonzero_degrees) > 10:
                # Log-log fit
                log_degrees = np.log(nonzero_degrees)
                # Not perfect but gives rough estimate
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    np.log(np.arange(1, len(nonzero_degrees) + 1)),
                    log_degrees,
                )
                metrics["degree_power_law_exponent"] = -slope
                metrics["degree_power_law_r2"] = r_value ** 2
        except Exception as e:
            logger.debug(f"Power law fit failed: {e}")

        logger.info("Network metrics computed successfully")
        return metrics

    def detect_communities(self, algorithm: str = "louvain") -> Dict[int, List[str]]:
        """
        Detect communities (clusters of densely-connected proteins).

        Args:
            algorithm: Community detection method:
                - "louvain": Fast modularity optimization (default, recommended)
                - "girvan_newman": Hierarchical edge-betweenness (slow, exact)
                - "label_propagation": Fast approximate (non-deterministic)

        Returns:
            Dict mapping community_id → list of protein node names
        """
        logger.info(f"Detecting communities using {algorithm}...")

        if algorithm == "louvain":
            # Louvain method (requires python-louvain package)
            partition = community_louvain.best_partition(
                self.G,
                weight="weight",
                random_state=self.random_seed,
            )
            # Convert: node -> community_id
            communities: Dict[int, List[str]] = {}
            for node, comm_id in partition.items():
                communities.setdefault(comm_id, []).append(node)

        elif algorithm == "girvan_newman":
            # Girvan-Newman (betweenness-based edge removal)
            comp = nx.algorithms.community.girvan_newman(self.G)
            # Take first partition (2 communities)
            top_level = next(comp)
            communities = {i: list(c) for i, c in enumerate(top_level)}

        elif algorithm == "label_propagation":
            communities_raw = nx.algorithms.community.label_propagation_communities(self.G)
            communities = {i: list(c) for i, c in enumerate(communities_raw)}

        else:
            raise ValueError(f"Unknown community algorithm: {algorithm}")

        logger.info(f"Detected {len(communities)} communities")
        # Log size distribution
        sizes = [len(v) for v in communities.values()]
        logger.debug(
            f"Community sizes: min={min(sizes)}, max={max(sizes)}, "
            f"mean={np.mean(sizes):.1f}, median={np.median(sizes)}"
        )

        return communities

    def find_hubs(
        self,
        top_n: int = 20,
        method: str = "degree",  # or "betweenness", "eigenvector"
    ) -> pd.DataFrame:
        """
        Identify hub proteins (key connectors in the network).

        Args:
            top_n: Number of top hubs to return
            method: Centrality measure to rank by

        Returns:
            DataFrame with columns: protein_id, centrality_score, degree, community_id
        """
        logger.info(f"Finding top {top_n} hubs by {method} centrality...")

        # Choose centrality
        if method == "degree":
            centrality = nx.degree_centrality(self.G)
        elif method == "betweenness":
            centrality = nx.betweenness_centrality(self.G, k=None, seed=self.random_seed)
        elif method == "eigenvector":
            try:
                centrality = nx.eigenvector_centrality_numpy(self.G, weight="weight")
            except Exception:
                logger.warning("Eigenvector failed, falling back to degree")
                centrality = nx.degree_centrality(self.G)
        else:
            raise ValueError(f"Unknown hub method: {method}")

        # Get communities for context
        communities = self.detect_communities(algorithm="louvain")
        node_to_comm = {}
        for comm_id, nodes in communities.items():
            for node in nodes:
                node_to_comm[node] = comm_id

        # Build hub table
        hub_data = []
        for node, score in sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:top_n]:
            hub_data.append({
                "protein_id": node,
                "centrality_score": round(score, 6),
                "degree": self.G.degree(node),
                "community_id": node_to_comm.get(node, -1),
            })

        df = pd.DataFrame(hub_data)
        logger.info(f"\nTop 10 hubs:\n{df.head(10).to_string()}")
        return df

    def compute_clustering(self) -> Dict[str, float]:
        """
        Compute clustering metrics (local and global).

        Returns:
            Dict with clustering coefficients
        """
        logger.info("Computing clustering coefficients...")

        results = {
            "global_transitivity": nx.transitivity(self.G),
            "avg_clustering": nx.average_clustering(self.G),
            "square_clustering": nx.square_clustering(self.G),  # dict per node
        }
        return results

    def get_network_summary(self) -> str:
        """
        Return human-readable summary of network properties.
        """
        metrics = self.compute_metrics()
        summary_lines = [
            "=" * 60,
            "NETWORK SUMMARY",
            "=" * 60,
            f"Nodes: {metrics['num_nodes']:,}",
            f"Edges: {metrics['num_edges']:,}",
            f"Density: {metrics['density']:.6f}",
            f"Connected components: {metrics.get('num_connected_components', 'N/A')}",
            f"Largest component: {metrics['largest_component_size']:,} "
            f"({metrics['largest_component_ratio']*100:.1f}%)",
            f"Diameter: {metrics.get('diameter', 'N/A')}",
            f"Avg shortest path: {metrics.get('avg_shortest_path_length', 'N/A'):.3f}"
            if metrics.get('avg_shortest_path_length') else "Avg shortest path: N/A",
            f"Clustering coefficient (global): {metrics['transitivity']:.4f}",
            f"Avg clustering: {metrics['avg_clustering_coefficient']:.4f}",
            f"Degree — mean: {metrics['degree_mean']:.2f}, max: {metrics['degree_max']}",
            f"Degree distribution exponent: {metrics.get('degree_power_law_exponent', 'N/A'):.2f}"
            if metrics.get('degree_power_law_exponent') else "Degree distribution: N/A",
        ]
        return "\n".join(summary_lines)

    def save_results(
        self,
        metrics: Dict,
        hubs: pd.DataFrame,
        communities: Dict,
        output_dir: Path,
        prefix: str = "network_analysis",
    ) -> None:
        """
        Save analysis results to disk.

        Args:
            metrics: Dict from compute_metrics()
            hubs: DataFrame from find_hubs()
            communities: Dict from detect_communities()
            output_dir: Output directory
            prefix: Filename prefix
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save metrics as JSON
        metrics_path = output_dir / f"{prefix}_metrics.json"
        with open(metrics_path, "w") as f:
            # Convert non-serializable values
            serializable = {}
            for k, v in metrics.items():
                if isinstance(v, (dict, list, str, int, float, bool, type(None))):
                    serializable[k] = v
                else:
                    serializable[k] = str(v)
            json.dump(serializable, f, indent=2)
        logger.info(f"Metrics saved → {metrics_path}")

        # Save hubs as CSV
        hubs_path = output_dir / f"{prefix}_hubs.csv"
        hubs.to_csv(hubs_path, index=False)
        logger.info(f"Hubs saved → {hubs_path}")

        # Save communities as JSON
        comm_path = output_dir / f"{prefix}_communities.json"
        with open(comm_path, "w") as f:
            # Convert to a more compact format: list of node lists
            comm_list = [nodes for nodes in communities.values()]
            json.dump(comm_list, f, indent=2)
        logger.info(f"Communities saved → {comm_path}")

        # Save summary text
        summary = self.get_network_summary()
        summary_path = output_dir / f"{prefix}_summary.txt"
        with open(summary_path, "w") as f:
            f.write(summary)
        logger.info(f"Summary saved → {summary_path}")


def load_network_from_csv(
    csv_path: Path,
    source_col: str = "protein1",
    target_col: str = "protein2",
    weight_col: Optional[str] = "combined_score",
    directed: bool = False,
) -> nx.Graph:
    """
    Load protein interaction network from a CSV/TSV file.

    File should have columns with protein identifiers and optionally a weight/score.

    Args:
        csv_path: Path to interaction table (CSV or TSV)
        source_col: Column name for protein A
        target_col: Column name for protein B
        weight_col: Column name for edge weight (confidence score). If None, unweighted.
        directed: Build directed graph? (Usually False for PPIs)

    Returns:
        NetworkX Graph

    Example:
        G = load_network_from_csv("data/processed/human_ppi.csv")
    """
    logger.info(f"Loading network from {csv_path}")

    df = pd.read_csv(csv_path, sep="\t" if csv_path.suffix == ".tsv" else ",")

    source_col = source_col or "protein1"
    target_col = target_col or "protein2"

    if source_col not in df.columns or target_col not in df.columns:
        raise ValueError(
            f"Required columns missing. Expected '{source_col}' and '{target_col}', "
            f"got: {list(df.columns)}"
        )

    G = nx.from_pandas_edgelist(
        df,
        source=source_col,
        target=target_col,
        edge_attr=weight_col if weight_col in df.columns else None,
        create_using=nx.DiGraph() if directed else nx.Graph(),
    )

    logger.info(
        f"Loaded graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges"
    )
    return G


if __name__ == "__main__":
    # Quick demo with a small network
    import sys

    from src.utils.logging import get_logger

    logger = get_logger(__name__)

    # Create a sample network
    G = nx.karate_club_graph()  # Classic test graph
    analyzer = GraphAnalyzer(G)
    metrics = analyzer.compute_metrics()
    hubs = analyzer.find_hubs(top_n=5)
    communities = analyzer.detect_communities()

    print(analyzer.get_network_summary())
    print("\nTop hubs:")
    print(hubs.to_string())
    print(f"\nCommunities: {len(communities)}")
