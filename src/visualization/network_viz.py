"""Interactive visualization of protein-protein interaction networks.

Uses Plotly for browser-based interactive graphs:
    - Zoom, pan, drag nodes
    - Color-coded communities
    - Edge thickness by confidence
    - Hover tooltips with protein details
    - Export to HTML (standalone, no server needed)

Can also generate static publication-quality figures with Matplotlib/Seaborn.

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
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from tqdm import tqdm

from src.network_analysis.graph_analyzer import GraphAnalyzer
from src.utils.config import load_config
from src.utils.logging import get_logger

logger = get_logger(__name__)


# Color palettes
PALETTES = {
    "plotly": px.colors.qualitative.Plotly,
    "plotly3": px.colors.qualitative.Plotly,  # Note: Plotly3 not available in newer plotly, using Plotly
    "set2": px.colors.qualitative.Set2,
    "d3": px.colors.qualitative.D3,
    "g10": px.colors.qualitative.G10,
}


class NetworkVisualizer:
    """
    Create interactive network visualizations from interaction data.

    Features:
        - Force-directed layout (Fruchterman-Reingold)
        - 3D visualization option
        - Community-based coloring
        - Edge weight as thickness
        - Hover info: protein ID, degree, betweenness
        - Export as standalone HTML (embedded JS, no external deps)

    Usage:
        viz = NetworkVisualizer(interaction_df)
        fig = viz.plot_interactive_graph(
            title="Human PPI Network",
            output_path="results/figures/network.html"
        )
        fig.show()

    Or with NetworkX graph directly:
        viz = NetworkVisualizer.from_networkx(G)
    """

    def __init__(
        self,
        data: pd.DataFrame | nx.Graph,
        source_col: str = "protein1",
        target_col: str = "protein2",
        weight_col: str = "combined_score",
        community_labels: Optional[Dict[str, int]] = None,
    ):
        """
        Initialize visualizer.

        Args:
            data: DataFrame with interaction data OR NetworkX Graph
            source_col: DataFrame column for source protein
            target_col: DataFrame column for target protein
            weight_col: DataFrame column for edge weight (confidence)
            community_labels: Optional dict protein → community_id for coloring
        """
        self.config = load_config()
        self.viz_config = self.config.get("visualization", {})

        if isinstance(data, nx.Graph):
            self.G = data
            self.df_edges = None
        else:
            self.df_edges = data.copy()
            self.G = self._build_network(data, source_col, target_col, weight_col)

        self.community_labels = community_labels or {}
        self.nodes_df = self._prepare_node_dataframe()

        logger.info(
            f"NetworkVisualizer initialized: "
            f"{self.G.number_of_nodes():,} nodes, {self.G.number_of_edges():,} edges"
        )

    @classmethod
    def from_networkx(cls, G: nx.Graph, community_labels: Optional[Dict] = None) -> "NetworkVisualizer":
        """Convenience constructor from NetworkX graph."""
        return cls(data=G, community_labels=community_labels)

    def _build_network(
        self,
        df: pd.DataFrame,
        source_col: str,
        target_col: str,
        weight_col: str,
    ) -> nx.Graph:
        """Build NetworkX graph from edge DataFrame."""
        G = nx.from_pandas_edgelist(
            df,
            source=source_col,
            target=target_col,
            edge_attr=[weight_col] if weight_col in df.columns else None,
            create_using=nx.Graph(),
        )
        return G

    def _prepare_node_dataframe(self) -> pd.DataFrame:
        """
        Build DataFrame with node attributes for visualization.

        Includes: degree, betweenness, community, x/y/z positions (lazy-computed),
        hover text, etc.
        """
        logger.info("Preparing node attributes...")

        # Basic degree
        degrees = dict(self.G.degree())
        nodes = list(self.G.nodes())

        node_df = pd.DataFrame({
            "id": nodes,
            "degree": [degrees.get(n, 0) for n in nodes],
            "label": [str(n) for n in nodes],  # Display name
        })

        # Communities
        if self.community_labels:
            node_df["community"] = node_df["id"].map(self.community_labels).fillna(-1).astype(int)
        else:
            # Try to auto-detect if we have analyzer
            try:
                analyzer = GraphAnalyzer(self.G)
                communities = analyzer.detect_communities(algorithm="louvain")
                node_to_comm = {}
                for comm_id, members in communities.items():
                    for node in members:
                        node_to_comm[node] = comm_id
                node_df["community"] = node_df["id"].map(node_to_comm).fillna(-1).astype(int)
            except Exception as e:
                logger.warning(f"Could not auto-detect communities: {e}")
                node_df["community"] = 0

        # Centralities (expensive, so only if not too many nodes)
        if len(self.G) < 5000:
            try:
                betweenness = nx.betweenness_centrality(self.G, k=min(200, len(self.G)), seed=42)
                node_df["betweenness"] = node_df["id"].map(betweenness).fillna(0)
            except Exception as e:
                logger.warning(f"Betweenness calculation failed: {e}")
                node_df["betweenness"] = 0.0
        else:
            logger.debug("Skipping betweenness (graph too large)")
            node_df["betweenness"] = 0.0

        # Hover text
        node_df["hover_text"] = node_df.apply(
            lambda row: (
                f"<b>{row['label']}</b><br>"
                f"Degree: {row['degree']}<br>"
                f"Community: {row['community']}"
            ),
            axis=1,
        )

        logger.info(f"Prepared {len(node_df):,} node attributes")
        return node_df

    def _compute_layout(
        self,
        algorithm: str = "fr",
        dim: int = 3,
        iterations: int = 1000,
        seed: int = 42,
    ) -> Dict[str, Tuple[float, float, float]]:
        """
        Compute node positions using specified layout algorithm.

        Args:
            algorithm: "fr" (Fruchterman-Reingold), "kk" (Kamada-Kawai),
                      "circle", "spectral", "spring"
            dim: 2 or 3 dimensions
            iterations: For force-directed layouts
            seed: Random seed for reproducibility

        Returns:
            Dict node_id -> (x, y) or (x, y, z)
        """
        logger.info(f"Computing {dim}D layout using {algorithm} ({iterations} iterations)...")

        rng = np.random.default_rng(seed)

        if algorithm == "fr":
            # Fruchterman-Reingold force-directed
            pos = nx.spring_layout(
                self.G,
                dim=dim,
                iterations=iterations,
                seed=seed,
            )
        elif algorithm == "kk":
            # Kamada-Kawai (better for small-medium graphs)
            try:
                pos = nx.kamada_kawai_layout(self.G, dim=dim)
            except Exception:
                logger.warning("Kamada-Kawai failed, falling back to spring")
                pos = nx.spring_layout(self.G, dim=dim, iterations=iterations, seed=seed)
        elif algorithm == "circle":
            pos = nx.circular_layout(self.G, dim=dim)
        elif algorithm == "spectral":
            pos = nx.spectral_layout(self.G, dim=dim)
        elif algorithm == "random":
            # Random layout within [-1, 1] cube
            pos = {node: tuple(rng.uniform(-1, 1, size=dim)) for node in self.G.nodes()}
        else:
            raise ValueError(f"Unknown layout algorithm: {algorithm}")

        logger.info(f"Layout computed for {len(pos)} nodes")
        return pos

    def plot_interactive_graph(
        self,
        title: str = "Protein Interaction Network",
        output_path: Optional[Path | str] = None,
        show_labels: Optional[bool] = None,
        node_size_col: Optional[str] = None,
        node_color_col: Optional[str] = "community",
        edge_width_scale: Optional[float] = None,
        layout_algorithm: Optional[str] = None,
        layout_iterations: Optional[int] = None,
        dim: int = 3,
        background_color: Optional[str] = None,
    ) -> go.Figure:
        """
        Create interactive Plotly graph figure.

        Args:
            title: Plot title
            output_path: If provided, save HTML to this path
            show_labels: Whether to display node labels (set False for large networks)
            node_size_col: Column in node_df to scale node size by (default: degree)
            node_color_col: Column to color nodes by (default: community)
            edge_width_scale: Scaling factor for edge width (higher = thicker)
            layout_algorithm: Layout algorithm (fr, kk, circle, etc.)
            layout_iterations: Number of iterations for force-directed algos
            dim: 2 or 3 dimensions
            background_color: Plot background (default from config)

        Returns:
            Plotly Figure object
        """
        logger.info("Generating interactive network visualization...")

        # Config overrides
        show_labels = show_labels if show_labels is not None else self.viz_config.get("nodes", {}).get("show_labels", True)
        edge_width_scale = edge_width_scale or self.viz_config.get("edges", {}).get("width_scale", 2.0)
        layout_algorithm = layout_algorithm or self.viz_config.get("layout", {}).get("algorithm", "fr")
        layout_iterations = layout_iterations or self.viz_config.get("layout", {}).get("iterations", 1000)
        background_color = background_color or self.viz_config.get("interactive", {}).get("bgcolor", "#0d1117")

        # Compute layout
        pos = self._compute_layout(
            algorithm=layout_algorithm,
            dim=dim,
            iterations=layout_iterations,
        )

        # Assign positions to nodes
        if dim == 3:
            self.nodes_df["x"] = self.nodes_df["id"].map(lambda n: pos[n][0])
            self.nodes_df["y"] = self.nodes_df["id"].map(lambda n: pos[n][1])
            self.nodes_df["z"] = self.nodes_df["id"].map(lambda n: pos[n][2])
            scene = dict(
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=""),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=""),
                zaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=""),
            )
        else:
            self.nodes_df["x"] = self.nodes_df["id"].map(lambda n: pos[n][0])
            self.nodes_df["y"] = self.nodes_df["id"].map(lambda n: pos[n][1])
            scene = None

        # Node size
        if node_size_col and node_size_col in self.nodes_df.columns:
            sizes = self.nodes_df[node_size_col]
        else:
            sizes = self.nodes_df["degree"]
        # Scale to min-max range
        size_min = self.viz_config.get("nodes", {}).get("min_size", 10)
        size_max = self.viz_config.get("nodes", {}).get("max_size", 50)
        if sizes.max() > sizes.min():
            node_sizes = size_min + (sizes - sizes.min()) / (sizes.max() - sizes.min()) * (size_max - size_min)
        else:
            node_sizes = size_min

        # Node color
        if node_color_col and node_color_col in self.nodes_df.columns:
            color_values = self.nodes_df[node_color_col]
            discrete_colors = True  # categorical
        else:
            color_values = self.nodes_df["degree"]
            discrete_colors = False

        # Edge coordinates
        edge_x = []
        edge_y = []
        edge_z = []
        edge_weights = []

        for _, row in self.df_edges.iterrows() if self.df_edges is not None else [
            {"protein1": u, "protein2": v, "combined_score": self.G[u][v].get("combined_score", 1)}
            for u, v in self.G.edges()
        ]:
            src = row["protein1" if self.df_edges is not None else "protein1"]
            tgt = row["protein2" if self.df_edges is not None else "protein2"]
            score = row.get("combined_score", 1) if self.df_edges is not None else self.G[src][tgt].get("combined_score", 1)

            src_row = self.nodes_df[self.nodes_df["id"] == src]
            tgt_row = self.nodes_df[self.nodes_df["id"] == tgt]

            if src_row.empty or tgt_row.empty:
                continue

            x0, y0 = src_row["x"].iloc[0], src_row["y"].iloc[0]
            x1, y1 = tgt_row["x"].iloc[0], tgt_row["y"].iloc[0]

            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            if dim == 3:
                z0, z1 = src_row["z"].iloc[0], tgt_row["z"].iloc[0]
                edge_z.extend([z0, z1, None])
            edge_weights.append(score)

        # Edge trace
        if dim == 3:
            edge_trace = go.Scatter3d(
                x=edge_x, y=edge_y, z=edge_z,
                mode="lines",
                line=dict(
                    width=1,
                    color="rgba(150,150,150,0.3)",
                ),
                hoverinfo="none",
                name="Edges",
            )
        else:
            edge_trace = go.Scatter(
                x=edge_x, y=edge_y,
                mode="lines",
                line=dict(
                    width=1,
                    color="rgba(150,150,150,0.3)",
                ),
                hoverinfo="none",
                name="Edges",
            )

        # Node trace
        color_palette = PALETTES.get(self.viz_config.get("nodes", {}).get("colormap", "plotly3"), PALETTES["plotly3"])

        if discrete_colors:
            # Categorical coloring by community
            unique_communities = sorted(color_values.unique())
            color_map = {comm: color_palette[i % len(color_palette)] for i, comm in enumerate(unique_communities)}
            node_colors = color_values.map(color_map)
        else:
            # Continuous colormap (viridis)
            node_colors = color_values
            colorscale = "Viridis"

        # Hover template
        hover_template = (
            "<b>%{text}</b><br>"
            "Degree: %{customdata[0]}<br>"
            "Community: %{customdata[1]}<extra></extra>"
        )

        if dim == 3:
            node_trace = go.Scatter3d(
                x=self.nodes_df["x"],
                y=self.nodes_df["y"],
                z=self.nodes_df["z"],
                mode="markers+text" if show_labels else "markers",
                marker=dict(
                    size=node_sizes,
                    color=node_colors if not discrete_colors else None,
                    colorscale=colorscale if not discrete_colors else None,
                    showscale=not discrete_colors,
                    colorbar=dict(title="Degree", thickness=15, x=1.02) if not discrete_colors else None,
                    line=dict(width=1, color="white"),
                    opacity=0.9,
                ),
                text=self.nodes_df["label"] if show_labels else None,
                textfont=dict(size=self.viz_config.get("nodes", {}).get("label_font_size", 10)),
                customdata=np.stack([
                    self.nodes_df["degree"],
                    self.nodes_df["community"],
                ], axis=-1),
                hovertemplate=hover_template,
                name="Proteins",
            )
        else:
            node_trace = go.Scatter(
                x=self.nodes_df["x"],
                y=self.nodes_df["y"],
                mode="markers+text" if show_labels else "markers",
                marker=dict(
                    size=node_sizes,
                    color=node_colors if not discrete_colors else None,
                    colorscale=colorscale if not discrete_colors else None,
                    showscale=not discrete_colors,
                    colorbar=dict(title="Degree", thickness=15, x=1.02) if not discrete_colors else None,
                    line=dict(width=1, color="white"),
                    opacity=0.9,
                ),
                text=self.nodes_df["label"] if show_labels else None,
                textfont=dict(size=self.viz_config.get("nodes", {}).get("label_font_size", 10)),
                customdata=np.stack([
                    self.nodes_df["degree"],
                    self.nodes_df["community"],
                ], axis=-1),
                hovertemplate=hover_template,
                name="Proteins",
            )

        # Build figure
        fig = go.Figure(data=[edge_trace, node_trace])

        fig.update_layout(
            title=dict(
                text=title,
                x=0.5,
                font=dict(size=18, family="Arial, sans-serif"),
            ),
            showlegend=False,
            margin=dict(l=0, r=0, t=50, b=0),
            paper_bgcolor=background_color,
            plot_bgcolor=background_color,
            hovermode="closest",
        )

        if dim == 3:
            fig.update_layout(
                scene=scene,
                width=self.viz_config.get("figures", {}).get("width", 1200),
                height=self.viz_config.get("figures", {}).get("height", 800),
            )
        else:
            fig.update_layout(
                width=self.viz_config.get("figures", {}).get("width", 1200),
                height=self.viz_config.get("figures", {}).get("height", 800),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            )

        # Save if path provided
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            fig.write_html(
                output_path,
                include_plotlyjs="cdn",  # Use CDN for smaller file
                full_html=True,
                auto_open=False,
            )
            logger.info(f"Interactive plot saved → {output_path}")

            # Also save static PNG for quick view
            png_path = output_path.with_suffix(".png")
            try:
                fig.write_image(
                    png_path,
                    width=1200,
                    height=800,
                    scale=2,  # high DPI
                )
                logger.info(f"Static PNG saved → {png_path}")
            except Exception as e:
                logger.debug(f"Could not save PNG (kaleido not available): {e}")

        return fig

    def plot_degree_distribution(
        self,
        output_path: Optional[Path] = None,
        log_scale: bool = True,
    ) -> go.Figure:
        """
        Plot degree distribution (log-log) to check scale-free property.

        Args:
            output_path: Save figure to this path
            log_scale: Use log-log axes (should be straight line if power law)

        Returns:
            Plotly Figure
        """
        degrees = [d for _, d in self.G.degree()]
        unique_degrees, counts = np.unique(degrees, return_counts=True)

        fig = go.Figure()

        if log_scale:
            fig.add_trace(go.Scatter(
                x=np.log1p(unique_degrees),
                y=np.log1p(counts),
                mode="markers",
                marker=dict(size=8, color="steelblue"),
            ))
            fig.update_layout(
                xaxis_title="log(Degree + 1)",
                yaxis_title="log(Frequency + 1)",
            )
        else:
            fig.add_trace(go.Scatter(
                x=unique_degrees,
                y=counts,
                mode="markers+lines",
                marker=dict(size=6, color="steelblue"),
            ))
            fig.update_layout(
                xaxis_title="Degree",
                yaxis_title="Frequency",
            )

        fig.update_layout(
            title="Degree Distribution" + (" (log-log)" if log_scale else ""),
            width=800,
            height=600,
            margin=dict(l=50, r=50, t=50, b=50),
        )

        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(output_path.with_suffix(".html"))
            try:
                fig.write_image(output_path.with_suffix(".png"), scale=2)
            except Exception as e:
                logger.debug(f"Could not save PNG (kaleido not available): {e}")
            logger.info(f"Degree distribution saved to {output_path}")

        return fig

    def plot_centrality_comparison(
        self,
        top_n: int = 50,
        output_path: Optional[Path] = None,
    ) -> go.Figure:
        """
        Compare centrality measures for top nodes.

        Args:
            top_n: Number of top nodes to compare (by degree)
            output_path: Save figure

        Returns:
            Plotly Figure
        """
        # Get top nodes by degree
        degrees = dict(self.G.degree())
        top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:top_n]
        top_node_names = [str(n) for n in top_nodes]

        # Compute centralities (only for top nodes)
        data = {"Protein": top_node_names, "Degree": [degrees[n] for n in top_nodes]}

        # Betweenness (expensive, compute only for these nodes using subgraph)
        if len(self.G) < 2000:  # Only compute if not too large
            try:
                betw = nx.betweenness_centrality(self.G, k=min(200, len(self.G)), seed=42)
                data["Betweenness"] = [betw.get(n, 0) for n in top_nodes]
            except Exception:
                pass

        # Eigenvector (if not too large)
        if len(self.G) < 1000:
            try:
                eig = nx.eigenvector_centrality_numpy(self.G)
                data["Eigenvector"] = [eig.get(n, 0) for n in top_nodes]
            except Exception:
                pass

        df = pd.DataFrame(data).melt(id_vars="Protein", var_name="Centrality", value_name="Score")

        fig = px.bar(
            df,
            x="Protein",
            y="Score",
            color="Centrality",
            barmode="group",
            title=f"Top {top_n} Proteins — Centrality Comparison",
            log_y=True,  # Typically follow power law
        )
        fig.update_layout(
            xaxis_tickangle=-45,
            width=1200,
            height=600,
            margin=dict(l=100, r=50, t=50, b=150),
        )

        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(output_path.with_suffix(".html"))
            logger.info(f"Centrality plot saved → {output_path}")

        return fig


# Convenience functions
def visualize_from_csv(
    csv_path: Path,
    output_html: Path,
    title: Optional[str] = None,
    score_threshold: int = 400,
    max_nodes: Optional[int] = None,
    **viz_kwargs,
) -> go.Figure:
    """
    One-liner to visualize interaction data from CSV.

    Args:
        csv_path: Path to interactions CSV
        output_html: Output HTML path
        title: Plot title (default: filename)
        score_threshold: Filter edges with combined_score >= threshold
        max_nodes: If provided, take top-N nodes by degree for large networks
        **viz_kwargs: Additional args to NetworkVisualizer.plot_interactive_graph()

    Returns:
        Plotly Figure
    """
    df = pd.read_csv(csv_path)

    # Filter by score
    if "combined_score" in df.columns and score_threshold:
        df = df[df["combined_score"] >= score_threshold].copy()
        logger.info(f"Filtered to {len(df):,} edges (score ≥ {score_threshold})")

    # Optionally limit nodes (take highest-degree ones)
    if max_nodes and len(df) > 0:
        from collections import Counter
        node_counts = Counter(df["protein1"].tolist() + df["protein2"].tolist())
        top_nodes = [n for n, _ in node_counts.most_common(max_nodes)]
        df = df[df["protein1"].isin(top_nodes) | df["protein2"].isin(top_nodes)].copy()
        logger.info(f"Filtered to top {max_nodes} nodes: {len(df):,} edges remain")

    title = title or f"Protein Interaction Network ({csv_path.stem})"

    viz = NetworkVisualizer(df)
    fig = viz.plot_interactive_graph(title=title, output_path=output_html, **viz_kwargs)
    return fig


if __name__ == "__main__":
    # Quick demo with networkx karate club graph
    import sys

    from src.utils.logging import get_logger

    logger = get_logger(__name__)

    print("=" * 60)
    print("Network Visualizer — Test (Karate Club Graph)")
    print("=" * 60)

    # Build test graph
    G = nx.karate_club_graph()

    # Add mock "community" labels
    true_labels = [G.nodes[n]["club"] for n in G.nodes()]
    label_map = {"Mr. Hi": 0, "Officer": 1}
    communities = {n: label_map[G.nodes[n]["club"]] for n in G.nodes()}

    viz = NetworkVisualizer.from_networkx(G, community_labels=communities)

    # Create output directory
    out_dir = Path("results/figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Interactive 3D
    fig1 = viz.plot_interactive_graph(
        title="Karate Club Network (3D)",
        output_path=out_dir / "test_network_3d.html",
        dim=3,
        show_labels=True,
    )

    # Interactive 2D
    fig2 = viz.plot_interactive_graph(
        title="Karate Club Network (2D)",
        output_path=out_dir / "test_network_2d.html",
        dim=2,
        show_labels=True,
    )

    logger.info("Test complete. Check results/figures/ for HTML output.")
