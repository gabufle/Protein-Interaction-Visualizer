"""Prediction module: machine learning for PPI prediction.

Provides:
    - InteractionPredictor: full ML pipeline (train, evaluate, predict)
    - load_and_predict: load saved model and predict new pairs

Author: [Your Name]
"""

from .interaction_predictor import InteractionPredictor, load_and_predict

__all__ = ["InteractionPredictor", "load_and_predict"]
