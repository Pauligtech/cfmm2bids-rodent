import logging
import sys
from pathlib import Path


def extract_subject_session_from_path(file_path):
    """
    Extract subject and session IDs from a file path.

    Expects path components like 'sub-{subject}' and 'ses-{session}' in directory names.
    Ignores the filename to avoid extracting from filenames that contain sub-/ses- patterns.

    Args:
        file_path: Path to the file

    Returns:
        tuple: (subject, session) or (None, None) if not found
    """
    path_obj = Path(file_path)
    # Use only directory parts, not the filename
    parts = path_obj.parent.parts if path_obj.parent else []
    subject = None
    session = None
    for part in parts:
        if part.startswith("sub-") and not subject:
            # Only extract if we haven't found one yet (take the first occurrence)
            subject = part.replace("sub-", "")
        elif part.startswith("ses-") and not session:
            # Only extract if we haven't found one yet (take the first occurrence)
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

    # Attach handler to THIS logger
    logger.addHandler(handler)

    # ⭐ IMPORTANT: attach handler to the ROOT logger too
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.addHandler(handler)

    if log_file:
        # Redirect stdout and stderr to the same log file
        sys.stdout = open(log_file, "a")
        sys.stderr = open(log_file, "a")

    return logger
