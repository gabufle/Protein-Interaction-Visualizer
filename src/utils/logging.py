import logging
import os
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path

from colorlog import ColoredFormatter

LOG_DIR = Path.home() / ".hermes" / "logs"
LOG_FILE = LOG_DIR / "protein_vis.log"

# Get project root
PROJECT_ROOT = Path(__file__).resolve().parents[2]

# Environment variable can override log level
DEFAULT_LOG_LEVEL = os.environ.get("PROTEIN_VIS_LOG_LEVEL", "INFO").upper()

# Log format (plain for file, colored for console)
FILE_FORMAT = logging.Formatter(
    fmt="%(asctime)s [%(levelname)-8s] %(name)s:%(lineno)d — %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
COLOR_FORMAT = ColoredFormatter(
    fmt="%(log_color)s%(asctime)s [%(levelname)-8s]%(reset)s %(name)s:%(lineno)d — %(message)s",
    datefmt="%H:%M:%S",
    log_colors={
        "DEBUG": "cyan",
        "INFO": "green",
        "WARNING": "yellow",
        "ERROR": "red",
        "CRITICAL": "bold_red",
    },
)


def _ensure_log_dir():
    """Make sure log directory exists."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)


def get_logger(name: str) -> logging.Logger:
    """
    Return a configured logger with both file and console handlers.

    Usage:
        logger = get_logger(__name__)
        logger.info("Hello world")
    """
    logger = logging.getLogger(name)
    logger.setLevel(DEFAULT_LOG_LEVEL)

    # Avoid adding handlers multiple times
    if logger.hasHandlers():
        return logger

    _ensure_log_dir()

    # File handler — full details, rotating
    file_handler = RotatingFileHandler(
        LOG_FILE,
        maxBytes=5 * 1024 * 1024,  # 5 MB
        backupCount=3,
        encoding='utf-8',
    )
    file_handler.setLevel(logging.DEBUG)  # File always gets full log
    file_handler.setFormatter(FILE_FORMAT)
    logger.addHandler(file_handler)

    # Console handler — colored output
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(DEFAULT_LOG_LEVEL)
    console_handler.setFormatter(COLOR_FORMAT)
    logger.addHandler(console_handler)

    return logger


def _get_effective_log_level() -> str:
    """Determine log level from env var or default."""
    env_level = os.environ.get("PROTEIN_VIS_LOG_LEVEL", "").upper()
    return env_level if env_level in ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL") else "INFO"


def silence_third_party_loggers(level: str = "WARNING") -> None:
    """
    Silence noisy loggers from third-party libraries.

    Call at startup before importing heavy libs.
    """
    noisy_loggers = [
        "matplotlib",
        "matplotlib.font_manager",
        "PIL",
        "plotly",
        "urllib3",
        "requests",
        "tqdm",
        "joblib",
        "sklearn",
        "community",
        "networkx",
    ]
    for name in noisy_loggers:
        logging.getLogger(name).setLevel(level)

    # Also silence warnings module (not a logger but prints to stderr)
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)


def set_global_log_level(level: str) -> None:
    """
    Change log level for all handlers on the fly.

    Args:
        level: DEBUG, INFO, WARNING, ERROR, CRITICAL
    """
    level = level.upper()
    os.environ["PROTEIN_VIS_LOG_LEVEL"] = level
    logger = logging.getLogger()
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setLevel(level)


# Convenience: call once at startup to set sensible defaults
def init_logging(level: str = None) -> logging.Logger:
    """
    Initialize logging for the whole project.

    Returns the root logger.
    """
    level = level or _get_effective_log_level()
    _ensure_log_dir()

    # Configure root logger
    root = logging.getLogger()
    root.setLevel(level)
    root.handlers.clear()

    # File handler
    fh = RotatingFileHandler(LOG_FILE, maxBytes=5*1024*1024, backupCount=3, encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(FILE_FORMAT)
    root.addHandler(fh)

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(COLOR_FORMAT)
    root.addHandler(ch)

    # Quiet third parties
    silence_third_party_loggers()

    root.info("=" * 60)
    root.info("Protein Interaction Visualizer — logging initialized")
    root.info(f"Log level: {level}")
    root.info(f"Log file: {LOG_FILE}")
    root.info("=" * 60)

    return root
