"""Configuration management for the protein interaction project.

Handles loading, validation, and access to project configuration.
Supports merging config.yaml with environment variables for secrets.

Usage:
    from src.utils.config import load_config
    config = load_config()
    score_threshold = config.data.string.score_threshold

Author: [Your Name]
Date: 2025-05-14
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional
import yaml
from dotenv import load_dotenv


# Project root: assume this file is at src/utils/config.py
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def load_config(
    config_path: Optional[Path] = None,
    env_path: Optional[Path] = None,
    overrides: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Load and merge configuration from multiple sources.

    Priority order (high to low):
        1. Function overrides parameter
        2. Environment variables (e.g. STRING_API_KEY)
        3. .env file (local secrets)
        4. config.yaml (default configuration)
        5. config_local.yaml (local overrides, gitignored)

    Args:
        config_path: Path to main config.yaml file (default: project_root/config.yaml)
        env_path: Path to .env file (default: project_root/.env)
        overrides: Dict of config values to override (useful for scripts)

    Returns:
        Dict containing merged configuration

    Raises:
        FileNotFoundError: If config.yaml is missing
        yaml.YAMLError: If config file is malformed
    """
    # Defaults
    if config_path is None:
        config_path = PROJECT_ROOT / "config.yaml"
    if env_path is None:
        env_path = PROJECT_ROOT / ".env"

    # 1. Load .env file (sets environment variables)
    if env_path.exists():
        load_dotenv(env_path, override=True)
    else:
        # .env is optional; warn if expected keys are missing later
        pass

    # 2. Load base YAML config
    if not config_path.exists():
        raise FileNotFoundError(
            f"Configuration file not found at {config_path}. "
            "Copy config.example.yaml to config.yaml and edit as needed."
        )

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # 3. Load local overrides (gitignored, for personal tweaks)
    local_config_path = PROJECT_ROOT / "config_local.yaml"
    if local_config_path.exists():
        with open(local_config_path, "r") as f:
            local_overrides = yaml.safe_load(f)
            _deep_update(config, local_overrides)

    # 4. Apply environment variable overrides for API keys
    # These are sensitive so they're top priority
    api_keys_section = config.get("api_keys", {})
    for key_name in api_keys_section:
        env_var = key_name.upper()  # STRING_API_KEY, BIOGRID_USERNAME, etc.
        if env_var in os.environ:
            api_keys_section[key_name] = os.environ[env_var]
            # Also expose at top level for convenience
            config["api_keys"][key_name] = os.environ[env_var]

    # 5. Apply function-level overrides
    if overrides:
        _deep_update(config, overrides)

    # 6. Validate required fields
    _validate_config(config)

    # Attach project_root for convenience
    config["_project_root"] = str(PROJECT_ROOT)

    return config


def _deep_update(base: Dict[str, Any], updates: Dict[str, Any]) -> None:
    """Recursively merge updates into base dict (in-place)."""
    for key, value in updates.items():
        if (
            key in base
            and isinstance(base[key], dict)
            and isinstance(value, dict)
        ):
            _deep_update(base[key], value)
        else:
            base[key] = value


def _validate_config(config: Dict[str, Any]) -> None:
    """Ensure required configuration values are present."""
    required_top_level = ["data", "network", "visualization", "prediction", "output"]
    for field in required_top_level:
        if field not in config:
            raise ValueError(f"Missing required config field: '{field}'")

    # Validate data section
    data_cfg = config["data"]
    if "string" not in data_cfg:
        raise ValueError("Config missing 'data.string' section")
    if "default_organism" not in data_cfg["string"]:
        raise ValueError("Config missing 'data.string.default_organism'")

    # Validate numeric thresholds
    score_thresh = data_cfg["string"].get("score_threshold", 400)
    if not (0 <= score_thresh <= 1000):
        raise ValueError(
            f"STRING score_threshold must be 0-1000, got {score_thresh}"
        )

    # Validate paths exist
    paths = config.get("paths", {})
    for path_key in ["data_raw", "data_processed", "results_figures"]:
        if path_key in paths:
            path_val = paths[path_key]
            full_path = PROJECT_ROOT / path_val
            full_path.mkdir(parents=True, exist_ok=True)


def get_log_level() -> str:
    """Get configured logging level from config or environment."""
    env_level = os.environ.get("PROTEIN_VIS_LOG_LEVEL", "").upper()
    if env_level in {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}:
        return env_level

    try:
        from src.utils.logging import get_logger
        config = load_config()
        return config.get("debug", {}).get("log_level", "INFO").upper()
    except Exception:
        return "INFO"


def get_random_seed() -> int:
    """Get random seed from config for reproducibility."""
    try:
        config = load_config()
        seed = config.get("debug", {}).get("random_seed", 42)
        if not isinstance(seed, int):
            raise ValueError(f"random_seed must be int, not {type(seed)}")
        return seed
    except Exception:
        return 42  # fallback


# Convenience: common config accessors
def get_organism_taxid() -> str:
    """Get default organism taxid (e.g., '9606' for human)."""
    config = load_config()
    return str(config["data"]["string"]["default_organism"])


def get_score_threshold() -> int:
    """Get STRING interaction score threshold (0-1000)."""
    config = load_config()
    return int(config["data"]["string"]["score_threshold"])


if __name__ == "__main__":
    # Test config loading
    import sys

    try:
        cfg = load_config()
        print("✅ Configuration loaded successfully")
        print(f"Project root: {cfg['_project_root']}")
        print(f"Organism: {get_organism_taxid()}")
        print(f"Score threshold: {get_score_threshold()}")
    except Exception as e:
        print(f"❌ Config error: {e}", file=sys.stderr)
        sys.exit(1)
