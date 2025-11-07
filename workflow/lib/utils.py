import logging
import sys
from pathlib import Path


def extract_subject_session_from_path(file_path):
    """
    Extract subject and session IDs from a file path.

    Expects path components like 'sub-{subject}' and 'ses-{session}'.

    Args:
        file_path: Path to the file

    Returns:
        tuple: (subject, session) or (None, None) if not found
    """
    parts = Path(file_path).parts
    subject = None
    session = None
    for part in parts:
        if part.startswith("sub-"):
            subject = part.replace("sub-", "")
        elif part.startswith("ses-"):
            session = part.replace("ses-", "")
    return subject, session


def setup_logger(log_file=None):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    if log_file:
        handler = logging.FileHandler(log_file)
    else:
        handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if log_file:
        # Redirect stdout and stderr to the same log file
        sys.stdout = open(log_file, "a")
        sys.stderr = open(log_file, "a")

    return logger
