"""Machine learning pipeline for predicting novel protein-protein interactions.

Uses network topology features + optional biological features to train
a binary classifier: will two proteins interact?

Approach:
    1. Positive examples: existing interactions from PPI database
    2. Negative examples: randomly sampled non-interacting protein pairs
    3. Features: topological (node2vec, common neighbors, Adamic-Adar, etc.)
                  biological (optional: sequence similarity, GO overlap)
    4. Model: scikit-learn RandomForest / GradientBoosting / Logistic Regression
    5. Evaluation: AUROC, Precision-Recall, confusion matrix
    6. Save model for reuse

Author: [Your Name]
Date: 2025-05-14
"""

import json
import pickle
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from joblib import dump, load
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    classification_report,
    confusion_matrix,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split
from tqdm import tqdm

from src.network_analysis.graph_analyzer import GraphAnalyzer
from src.utils.config import load_config, get_random_seed
from src.utils.logging import get_logger

logger = get_logger(__name__)


class InteractionPredictor:
    """
    Train ML models to predict protein-protein interactions.

    Pipeline:
        1. Prepare graph from interaction data
        2. Generate positive (real) and negative (random) protein pairs
        3. Extract features for each pair (topological + optional biological)
        4. Split into train/test
        5. Train classifier (Random Forest by default)
        6. Evaluate performance
        7. Save model for predictions on new protein pairs

    Usage:
        predictor = InteractionPredictor(G)
        X, y, pairs = predictor.prepare_features(neg_ratio=1.0)
        predictor.train(X, y)
        results = predictor.evaluate(X_test, y_test)
        predictor.save_models("models/")
    """

    # Available topological features
    TOPO_FEATURES = {
        "common_neighbors": "Number of shared neighbors",
        "jaccard_coefficient": "|N(u) ∩ N(v)| / |N(u) ∪ N(v)|",
        "adamic_adar": "Σ 1/log(|N(w)|) for shared neighbors w",
        "preferential_attachment": "|N(u)| × |N(v)|",
        "resource_allocation": "Σ 1/|N(w)| for shared neighbors w",
        "cn_soundarajan_hopcroft": "Common neighbors with community (if known)",
        "within_inter_cluster": "Connectivity within same cluster",
    }

    # Available biological features (require additional data)
    BIO_FEATURES = {
        "sequence_similarity": "BLAST/FASTA alignment score (requires sequences)",
        "go_jaccard": "Gene Ontology term overlap",
        "go_semantic_similarity": "GO term semantic similarity",
        "domain_overlap": "Pfam domain intersection",
        "coexpression": "Gene co-expression correlation (requires expression data)",
    }

    def __init__(
        self,
        graph: nx.Graph,
        communities: Optional[Dict] = None,
        protein_annotations: Optional[Dict[str, Dict]] = None,
        config: Optional[Dict] = None,
    ):
        """
        Initialize predictor.

        Args:
            graph: NetworkX graph of protein interactions (positive examples)
            communities: Optional dict node → community_id for community features
            protein_annotations: Optional dict protein → annotation dict
                (e.g., {"P00533": {"go_terms": [...], "domains": [...]}})
            config: Optional config override
        """
        self.G = graph
        self.communities = communities or {}
        self.annotations = protein_annotations or {}
        self.config = config or load_config()
        self.pred_config = self.config.get("prediction", {})
        self.random_seed = get_random_seed()

        self.model = None
        self.feature_names = []
        self.scaler = None  # For standardization
        self.training_results = {}

        logger.info(
            f"InteractionPredictor initialized on graph with "
            f"{self.G.number_of_nodes():,} nodes, {self.G.number_of_edges():,} edges"
        )

    def prepare_features(
        self,
        neg_ratio: float = 1.0,
        test_size: float = 0.2,
        val_size: float = 0.1,
        feature_set: Optional[List[str]] = None,
    ) -> Tuple[np.ndarray, np.ndarray, List[Tuple[str, str]]]:
        """
        Generate feature matrix and labels from graph.

        Args:
            neg_ratio: Number of negative samples per positive (1.0 = balanced)
            test_size: Fraction of data for test set
            val_size: Fraction of training data for validation
            feature_set: List of feature names to include. If None, uses
                        config-defined or all available topological features.

        Returns:
            X: Feature matrix (n_samples, n_features)
            y: Binary labels (1 = positive, 0 = negative)
            pairs: List of (protein_A, protein_B) tuples corresponding to rows

        Raises:
            ValueError: If graph is too small or has no edges
        """
        logger.info("Preparing features for ML...")

        # Determine feature set
        if feature_set is None:
            feature_set = list(self.TOPO_FEATURES.keys())
            # Optionally add biological features if annotations present
            if self.annotations and self.pred_config.get("features", {}).get("use_go_similarity", False):
                feature_set.append("go_jaccard")
            if self.annotations and self.pred_config.get("features", {}).get("use_domain_features", False):
                feature_set.append("domain_overlap")

        logger.info(f"Using features: {feature_set}")
        self.feature_names = feature_set

        # 1. Get positive pairs (existing edges)
        positives = list(self.G.edges())
        logger.info(f"Positive examples: {len(positives):,}")

        # 2. Generate negative pairs (non-edges)
        n_neg = int(len(positives) * neg_ratio)
        logger.info(f"Generating {n_neg:,} negative samples...")

        negatives = self._sample_negative_pairs(n_neg, exclude=positives)
        logger.info(f"Negative samples: {len(negatives):,}")

        # 3. Build feature vectors
        all_pairs = positives + negatives
        labels = np.array([1] * len(positives) + [0] * len(negatives))

        X = self._extract_features(all_pairs, feature_set)

        logger.info(
            f"Feature matrix: {X.shape[0]:,} samples × {X.shape[1]} features"
        )

        # 4. Train/val/test split
        # First split off test
        X_train, X_test, y_train, y_test, pairs_train, pairs_test = train_test_split(
            X, labels, all_pairs,
            test_size=test_size,
            random_state=self.random_seed,
            stratify=labels,
        )

        # Then split training into train/val
        if val_size > 0:
            val_ratio = val_size / (1 - test_size)
            X_train, X_val, y_train, y_val = train_test_split(
                X_train, y_train,
                test_size=val_ratio,
                random_state=self.random_seed,
                stratify=y_train,
            )
            logger.info(
                f"Split: train={len(X_train):,}, val={len(X_val):,}, test={len(X_test):,}"
            )
            self.val_data = (X_val, y_val)
        else:
            logger.info(f"Split: train={len(X_train):,}, test={len(X_test):,}")
            self.val_data = None

        self.train_data = (X_train, y_train)
        self.test_data = (X_test, y_test)
        self.pairs_test = pairs_test

        return X, labels, all_pairs

    def _sample_negative_pairs(
        self,
        n_samples: int,
        exclude: List[Tuple[str, str]],
    ) -> List[Tuple[str, str]]:
        """
        Randomly sample non-edges from the graph (protein pairs with no edge).

        Strategy: randomly sample node pairs and keep those not in graph.
        For large graphs, rejection sampling is efficient.

        Args:
            n_samples: Number of negative samples to generate
            exclude: Set of positive pairs to exclude

        Returns:
            List of (node1, node2) tuples
        """
        nodes = list(self.G.nodes())
        exclude_set = set(exclude)
        negatives = []

        attempts = 0
        max_attempts = n_samples * 10  # Prevent infinite loops on dense graphs

        with tqdm(total=n_samples, desc="Sampling negatives") as pbar:
            while len(negatives) < n_samples and attempts < max_attempts:
                u = np.random.choice(nodes)
                v = np.random.choice(nodes)

                if u == v:
                    continue  # no self-loops

                pair = (u, v) if u < v else (v, u)  # canonical ordering

                if pair in exclude_set:
                    continue

                if pair not in negatives:
                    negatives.append(pair)
                    pbar.update(1)

                attempts += 1

        if len(negatives) < n_samples:
            logger.warning(
                f"Could only generate {len(negatives):,}/{n_samples:,} negatives. "
                f"Graph may be too dense."
            )

        return negatives

    def _extract_features(
        self,
        pairs: List[Tuple[str, str]],
        feature_names: List[str],
    ) -> np.ndarray:
        """
        Compute feature vectors for protein pairs.

        Args:
            pairs: List of (node_A, node_B) tuples
            feature_names: Which features to compute

        Returns:
            Feature matrix (n_pairs, n_features)
        """
        logger.info(f"Extracting {len(feature_names)} features for {len(pairs):,} pairs...")

        features = []

        for u, v in tqdm(pairs, desc="Computing features"):
            vec = []

            # Topological features (using NetworkX)
            if "common_neighbors" in feature_names:
                cn = len(list(nx.common_neighbors(self.G, u, v)))
                vec.append(cn)

            if "jaccard_coefficient" in feature_names:
                try:
                    jc = list(nx.jaccard_coefficient(self.G, [(u, v)]))[0][2]
                    vec.append(jc)
                except:
                    vec.append(0.0)

            if "adamic_adar" in feature_names:
                try:
                    aa = list(nx.adamic_adar_index(self.G, [(u, v)]))[0][2]
                    vec.append(aa)
                except:
                    vec.append(0.0)

            if "preferential_attachment" in feature_names:
                pa = self.G.degree(u) * self.G.degree(v)
                vec.append(pa)

            if "resource_allocation" in feature_names:
                try:
                    ra = list(nx.resource_allocation_index(self.G, [(u, v)]))[0][2]
                    vec.append(ra)
                except:
                    vec.append(0.0)

            if "cn_soundarajan_hopcroft" in feature_names:
                if self.communities:
                    comm_u = self.communities.get(u, -1)
                    comm_v = self.communities.get(v, -1)
                    if comm_u == comm_v and comm_u != -1:
                        shared = []
                        for w in nx.common_neighbors(self.G, u, v):
                            if self.communities.get(w, -1) == comm_u:
                                shared.append(w)
                        vec.append(len(shared))
                    else:
                        vec.append(0)
                else:
                    vec.append(0)

            # Biological features (if available)
            if "go_jaccard" in feature_names:
                go_u = set(self.annotations.get(u, {}).get("go_terms", []))
                go_v = set(self.annotations.get(v, {}).get("go_terms", []))
                if go_u and go_v:
                    inter = len(go_u & go_v)
                    union = len(go_u | go_v)
                    vec.append(inter / union if union > 0 else 0.0)
                else:
                    vec.append(0.0)

            if "domain_overlap" in feature_names:
                dom_u = set(self.annotations.get(u, {}).get("domains", []))
                dom_v = set(self.annotations.get(v, {}).get("domains", []))
                if dom_u and dom_v:
                    inter = len(dom_u & dom_v)
                    union = len(dom_u | dom_v)
                    vec.append(inter / union if union > 0 else 0.0)
                else:
                    vec.append(0.0)

            features.append(vec)

        # Pad feature vectors to match requested feature count
        # If any features were skipped (e.g., community-based with no data),
        # ensure all requested features have a column by adding zeros as needed.
        expected_len = len(feature_names)
        for i, vec in enumerate(features):
            if len(vec) < expected_len:
                # Pad with zeros at the end (corresponding to skipped features)
                features[i].extend([0.0] * (expected_len - len(vec)))
            elif len(vec) > expected_len:
                # Truncate (shouldn't happen, but safe)
                features[i] = vec[:expected_len]

        return np.array(features, dtype=np.float32)

    def train(
        self,
        X_train: np.ndarray,
        y_train: np.ndarray,
        model_type: str = "random_forest",
        **model_kwargs,
    ) -> None:
        """
        Train the interaction predictor.

        Args:
            X_train: Training features
            y_train: Training labels
            model_type: "random_forest" (default), "gradient_boosting", "logistic_regression"
            **model_kwargs: Model-specific hyperparameters (overrides config)
        """
        logger.info(f"Training {model_type} model...")

        # Get model config
        model_cfg = self.pred_config.get("model", {})
        if model_type is None:
            model_type = model_cfg.get("algorithm", "random_forest")

        # Standardize features
        from sklearn.preprocessing import StandardScaler
        self.scaler = StandardScaler()
        X_train_scaled = self.scaler.fit_transform(X_train)

        # Initialize model
        random_state = model_kwargs.pop("random_state", self.random_seed)

        if model_type == "random_forest":
            self.model = RandomForestClassifier(
                n_estimators=model_kwargs.get("n_estimators", model_cfg.get("n_estimators", 200)),
                max_depth=model_kwargs.get("max_depth", model_cfg.get("max_depth", 20)),
                min_samples_split=model_kwargs.get("min_samples_split", model_cfg.get("min_samples_split", 10)),
                class_weight=model_kwargs.get("class_weight", model_cfg.get("class_weight", "balanced")),
                random_state=random_state,
                n_jobs=-1,  # Use all cores
                verbose=0,
            )
        elif model_type == "gradient_boosting":
            from sklearn.ensemble import GradientBoostingClassifier
            self.model = GradientBoostingClassifier(
                n_estimators=model_kwargs.get("n_estimators", 200),
                max_depth=model_kwargs.get("max_depth", 5),
                random_state=random_state,
            )
        elif model_type == "logistic_regression":
            self.model = LogisticRegression(
                penalty="l2",
                C=model_kwargs.get("C", 1.0),
                class_weight=model_kwargs.get("class_weight", "balanced"),
                random_state=random_state,
                max_iter=1000,
                n_jobs=-1,
            )
        else:
            raise ValueError(f"Unsupported model type: {model_type}")

        # Train
        self.model.fit(X_train_scaled, y_train)
        logger.info("Training complete")

        # Feature importance
        if hasattr(self.model, "feature_importances_"):
            importance_df = pd.DataFrame({
                "feature": self.feature_names,
                "importance": self.model.feature_importances_,
            }).sort_values("importance", ascending=False)
            logger.info("Feature importances:")
            for _, row in importance_df.iterrows():
                logger.info(f"  {row['feature']}: {row['importance']:.4f}")
            self.feature_importance = importance_df
        else:
            self.feature_importance = None

        self.model_type = model_type
        self.training_results["model_type"] = model_type
        self.training_results["n_features"] = X_train.shape[1]
        self.training_results["n_samples"] = X_train.shape[0]

    def evaluate(
        self,
        X_test: np.ndarray,
        y_test: np.ndarray,
        pairs_test: Optional[List[Tuple[str, str]]] = None,
    ) -> Dict:
        """
        Evaluate trained model on test set.

        Args:
            X_test: Test features
            y_test: Test labels
            pairs_test: Optional list of protein pairs for analysis

        Returns:
            Dict of metrics
        """
        if self.model is None:
            raise ValueError("Model not trained. Call .train() first.")

        logger.info("Evaluating model on test set...")

        # Scale features
        X_test_scaled = self.scaler.transform(X_test)

        # Predictions
        y_pred = self.model.predict(X_test_scaled)
        y_proba = self.model.predict_proba(X_test_scaled)[:, 1] if hasattr(self.model, "predict_proba") else y_pred

        # Metrics
        acc = accuracy_score(y_test, y_pred)
        roc_auc = roc_auc_score(y_test, y_proba)
        avg_prec = average_precision_score(y_test, y_proba)

        logger.info(f"Accuracy: {acc:.4f}")
        logger.info(f"ROC AUC: {roc_auc:.4f}")
        logger.info(f"Avg Precision: {avg_prec:.4f}")

        # Detailed report
        report = classification_report(y_test, y_pred, output_dict=True)
        cm = confusion_matrix(y_test, y_pred)

        results = {
            "accuracy": acc,
            "roc_auc": roc_auc,
            "average_precision": avg_prec,
            "classification_report": report,
            "confusion_matrix": cm.tolist(),
            "feature_importance": self.feature_importance.to_dict() if self.feature_importance is not None else None,
        }

        # Save to results
        self.training_results["test_metrics"] = {
            "accuracy": acc,
            "roc_auc": roc_auc,
            "average_precision": avg_prec,
        }

        # Optional: save predictions for error analysis
        if pairs_test:
            pred_df = pd.DataFrame({
                "protein_A": [p[0] for p in pairs_test],
                "protein_B": [p[1] for p in pairs_test],
                "true_label": y_test,
                "predicted": y_pred,
                "probability": y_proba,
            })
            results["predictions"] = pred_df

        return results

    def save_models(self, output_dir: Path, save_scaler: bool = True) -> None:
        """
        Save trained model and scaler to disk.

        Args:
            output_dir: Directory to save model files
            save_scaler: Save the StandardScaler for future feature scaling
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Save model
        model_path = output_dir / f"ppi_predictor_{self.model_type}_{timestamp}.joblib"
        dump(self.model, model_path)
        logger.info(f"Model saved → {model_path}")

        # Save scaler
        if save_scaler and self.scaler:
            scaler_path = output_dir / f"scaler_{timestamp}.joblib"
            dump(self.scaler, scaler_path)
            logger.info(f"Scaler saved → {scaler_path}")

        # Save metadata
        metadata = {
            "model_type": self.model_type,
            "feature_names": self.feature_names,
            "training_date": datetime.now().isoformat(),
            "n_features": self.model.n_features_in_ if hasattr(self.model, "n_features_in_") else len(self.feature_names),
            "model_path": str(model_path),
            "scaler_path": str(scaler_path) if save_scaler else None,
            "training_results": self.training_results,
            "config_snapshot": self.pred_config,
        }

        meta_path = output_dir / f"predictor_metadata_{timestamp}.json"
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=2, default=str)
        logger.info(f"Metadata saved → {meta_path}")

        # Also save latest symlink/copy for convenience
        latest_model = output_dir / "ppi_predictor_latest.joblib"
        latest_scaler = output_dir / "scaler_latest.joblib"
        import shutil
        shutil.copy(model_path, latest_model)
        if save_scaler:
            shutil.copy(scaler_path, latest_scaler)
        logger.info("Created 'latest' symlinks for easy loading")

    @classmethod
    def load_models(
        cls,
        model_path: Path,
        scaler_path: Optional[Path] = None,
        metadata_path: Optional[Path] = None,
        graph: Optional[nx.Graph] = None,
    ) -> "InteractionPredictor":
        """
        Load a saved predictor from disk.

        Args:
            model_path: Path to saved model (.joblib)
            scaler_path: Path to saved scaler (.joblib)
            metadata_path: Path to metadata JSON (optional, for feature names)
            graph: Graph to associate with (optional)

        Returns:
            InteractionPredictor instance with loaded model
        """
        model_path = Path(model_path)
        scaler_path = Path(scaler_path) if scaler_path else model_path.parent / "scaler_latest.joblib"

        model = load(model_path)
        scaler = load(scaler_path) if scaler_path.exists() else None

        # Load metadata if available
        if metadata_path is None:
            # Guess metadata path from model name
            meta_name = model_path.name.replace(".joblib", "_metadata.json")
            metadata_path = model_path.parent / meta_name

        feature_names = []
        if metadata_path and Path(metadata_path).exists():
            with open(metadata_path) as f:
                metadata = json.load(f)
            feature_names = metadata.get("feature_names", [])

        predictor = cls(graph or nx.Graph())
        predictor.model = model
        predictor.scaler = scaler
        predictor.feature_names = feature_names
        predictor.model_type = model_path.stem.split("_")[1] if "_" in model_path.stem else "unknown"

        logger.info(f"Loaded {predictor.model_type} model from {model_path}")
        return predictor

    def predict(
        self,
        pairs: List[Tuple[str, str]],
        threshold: float = 0.5,
    ) -> pd.DataFrame:
        """
        Predict interaction probability for new protein pairs.

        Args:
            pairs: List of (protein_A, protein_B) tuples
            threshold: Decision threshold (default 0.5)

        Returns:
            DataFrame with columns: protein_A, protein_B, probability, predicted
        """
        if self.model is None:
            raise ValueError("Model not loaded. Train a model or load from disk first.")

        logger.info(f"Predicting {len(pairs):,} protein pairs...")

        # Extract features (using the graph's current state)
        X = self._extract_features(pairs, self.feature_names)

        # Scale
        X_scaled = self.scaler.transform(X)

        # Predict
        proba = self.model.predict_proba(X_scaled)[:, 1]
        pred = (proba >= threshold).astype(int)

        results = pd.DataFrame({
            "protein_A": [p[0] for p in pairs],
            "protein_B": [p[1] for p in pairs],
            "probability": proba,
            "predicted": pred,
        })

        n_pos = (pred == 1).sum()
        logger.info(f"Predictions: {n_pos:,}/{len(pairs):,} predicted as interactions")

        return results

    def save_predictions(
        self,
        predictions: pd.DataFrame,
        output_path: Path,
        format: str = "csv",
    ) -> None:
        """Save predictions to disk."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == "csv":
            predictions.to_csv(output_path, index=False)
        elif format == "tsv":
            predictions.to_csv(output_path, sep="\t", index=False)
        elif format == "parquet":
            predictions.to_parquet(output_path, index=False)
        else:
            raise ValueError(f"Unsupported format: {format}")

        logger.info(f"Predictions saved → {output_path}")


# Utility: load trained predictor and apply to new data
def load_and_predict(
    model_dir: Path,
    pairs: List[Tuple[str, str]],
    graph: Optional[nx.Graph] = None,
    threshold: float = 0.5,
) -> pd.DataFrame:
    """
    One-liner: load saved model and predict on new pairs.

    Example:
        predictor = load_and_predict("models/", [("P00533", "P04637")])
        print(predictions)
    """
    model_path = Path(model_dir) / "ppi_predictor_latest.joblib"
    scaler_path = Path(model_dir) / "scaler_latest.joblib"
    meta_path = Path(model_dir) / "predictor_metadata_latest.json"

    predictor = InteractionPredictor.load_models(
        model_path=model_path,
        scaler_path=scaler_path,
        metadata_path=meta_path,
        graph=graph,
    )
    predictions = predictor.predict(pairs, threshold=threshold)
    return predictions


if __name__ == "__main__":
    import sys

    from src.utils.logging import get_logger

    logger = get_logger(__name__)

    print("=" * 60)
    print("Interaction Predictor — Test Mode")
    print("=" * 60)

    # Generate a small test graph
    G = nx.erdos_renyi_graph(n=200, p=0.05, seed=42)

    logger.info(f"Test graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    predictor = InteractionPredictor(G)

    # Prepare features
    X, y, pairs = predictor.prepare_features(neg_ratio=1.0, test_size=0.3)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42, stratify=y
    )

    # Train
    predictor.train(X_train, y_train, model_type="random_forest")

    # Evaluate
    results = predictor.evaluate(X_test, y_test)

    print("\nTest Results:")
    print(f"  Accuracy:  {results['accuracy']:.3f}")
    print(f"  ROC AUC:   {results['roc_auc']:.3f}")
    print(f"  Avg Prec:  {results['average_precision']:.3f}")

    # Save models
    out_dir = Path("models") / "test_run"
    predictor.save_models(out_dir)

    logger.info("✅ Test complete")
