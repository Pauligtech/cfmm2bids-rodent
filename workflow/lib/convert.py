"""Helper functions for DICOM to BIDS conversion with heudiconv."""

import ast
import logging
import shutil
import stat
import subprocess
from collections import defaultdict
from pathlib import Path
from pprint import pprint

import pandas as pd

logger = logging.getLogger(__name__)

root = logging.getLogger()
if not root.handlers:
    raise RuntimeError("No root handlers configured")

log_fh = root.handlers[0].stream


def find_tar_files(dicoms_dir: Path) -> list[Path]:
    """Find all tar files in the dicoms directory."""
    tar_files = list(dicoms_dir.glob("*.tar"))
    if not tar_files:
        tar_files = list(dicoms_dir.glob("*.tar.gz"))
    return sorted(tar_files)


def extract_study_uid_from_tar(tar_path: Path) -> str:
    """Extract the StudyInstanceUID from the tar filename or first DICOM."""
    # Try to get it from the tar filename first
    # cfmm2tar typically names files like: <StudyInstanceUID>.tar
    stem = tar_path.stem
    if stem.endswith(".tar"):
        stem = stem[:-4]
    return stem


def make_tree_writable(path: Path):
    """Recursively make files and dirs under `path` writable by owner."""
    for p in [path] + list(path.rglob("*")):
        try:
            mode = p.stat().st_mode
            p.chmod(mode | stat.S_IWUSR)
        except Exception as e:
            logger.info(f"⚠️ Could not make writable: {p} ({e})")


def run_heudiconv_for_study(
    tar_file: Path,
    temp_dir: Path,
    subject: str,
    session: str,
    heuristic: Path,
    dcmconfig_json: Path,
    heudiconv_options: str,
) -> dict[str, Path]:
    """
    Run heudiconv for a single tar file (study).

    Returns dict with paths to generated files.
    """
    study_uid = extract_study_uid_from_tar(tar_file)
    logger.info(f"Processing study {study_uid} from {tar_file.name}")

    # Create a temporary directory for this study
    study_temp_dir = Path(temp_dir) / f"study_{study_uid}"
    study_temp_dir.mkdir(parents=True, exist_ok=True)

    # Run heudiconv on this study's DICOMs
    bids_output = study_temp_dir / "bids"
    bids_output.mkdir(parents=True, exist_ok=True)

    cmd = [
        "heudiconv",
        "--files",
        str(tar_file),
        "-c",
        "dcm2niix",
        "-o",
        str(bids_output),
        "-ss",
        session,
        "-s",
        subject,
        "-f",
        str(heuristic),
        "--bids",
        "notop",
        "--dcmconfig",
        str(dcmconfig_json),
        "--overwrite",
    ]

    if heudiconv_options:
        cmd.extend(heudiconv_options.split())

    logger.info(f"Running heudiconv: {' '.join(cmd)}")
    result = subprocess.run(cmd, stdout=log_fh, stderr=log_fh, text=True)

    if result.returncode != 0:
        logger.error(f"heudiconv failed for study {study_uid}:")
        logger.error(result.stderr)
        raise RuntimeError(f"heudiconv failed for study {study_uid}")

    # Locate the heudiconv info files
    heudiconv_info_dir = (
        bids_output / ".heudiconv" / subject / f"ses-{session}" / "info"
    )

    auto_txt = heudiconv_info_dir / f"{subject}_ses-{session}.auto.txt"
    dicominfo_tsv = heudiconv_info_dir / f"dicominfo_ses-{session}.tsv"

    # Check that files exist
    for f in [auto_txt, dicominfo_tsv]:
        if not f.exists():
            raise FileNotFoundError(f"Expected heudiconv output file not found: {f}")

    return {
        "study_uid": study_uid,
        "bids_dir": bids_output / f"sub-{subject}" / f"ses-{session}",
        "auto_txt": auto_txt,
        "dicominfo_tsv": dicominfo_tsv,
    }


def merge_dicominfo_files(info_files: list[dict[str, Path]], output_tsv: Path) -> None:
    """
    Merge dicominfo.tsv files with series ID offsetting.

    Each study's series IDs are offset to avoid conflicts (only for numeric IDs).
    """
    logger.info("Merging dicominfo.tsv files...")

    all_dfs = []

    for _i, info in enumerate(info_files):
        df = pd.read_csv(info["dicominfo_tsv"], sep="\t")

        # Add study_uid column to track which study each series came from
        df["study_uid"] = info["study_uid"]
        all_dfs.append(df)

    merged_df = pd.concat(all_dfs, ignore_index=True)
    merged_df.to_csv(output_tsv, sep="\t", index=False)
    logger.info(f"Merged dicominfo saved to {output_tsv}")


def merge_auto_txt_files(info_files: list[dict[str, Path]], output_txt: Path) -> None:
    """
    Merge auto.txt files.

    Simply concatenate all auto.txt files with study separators.
    """
    logger.info("Merging auto.txt files...")

    merged = defaultdict(list)

    for _i, info in enumerate(info_files):
        with open(info["auto_txt"]) as f:
            text = f.read().strip()
        try:
            data = ast.literal_eval(text)
        except Exception as e:
            print(f"⚠️ Skipping info['auto_txt']: could not parse ({e})")
            continue

        for k, v in data.items():
            if not isinstance(v, list):
                continue
            merged[k].extend(v)

    # --- SAVE RESULT ---
    with open(output_txt, "w") as f:
        # pretty-print to make it readable (preserves valid Python syntax)
        pprint(dict(merged), stream=f, sort_dicts=False, width=120, compact=False)

    logger.info(f"Merged auto.txt saved to {output_txt}")


def merge_bids_directories(
    info_files: list[dict[str, Path]], output_bids_dir: Path
) -> None:
    """
    Merge BIDS directories from multiple studies.

    Copies all NIfTI and JSON files into a single output directory.
    """
    logger.info("Merging BIDS directories...")

    output_bids_dir.mkdir(parents=True, exist_ok=True)

    for info in info_files:
        src_dir = info["bids_dir"]
        if not src_dir.exists():
            logger.warning(f"BIDS directory not found: {src_dir}")
            continue

        # First make writable before copying, as perms preserved
        make_tree_writable(src_dir)

        # Copy all files from source to destination
        for src_file in src_dir.rglob("*"):
            if src_file.is_file():
                rel_path = src_file.relative_to(src_dir)
                dst_file = output_bids_dir / rel_path
                dst_file.parent.mkdir(parents=True, exist_ok=True)

                # clobber files

                shutil.copy2(src_file, dst_file)
                logger.debug(f"Copied {src_file} to {dst_file}")

    logger.info(f"Merged BIDS directory created at {output_bids_dir}")


def process_single_study_heudiconv(
    dicoms_dir: Path,
    subject: str,
    session: str,
    heuristic: Path,
    dcmconfig_json: Path,
    heudiconv_options: str,
    output_bids_dir: Path,
    output_info_dir: Path,
    temp_dir: Path,
) -> None:
    """
    Process a single tar file with heudiconv.

    Args:
        dicoms_dir: Directory containing tar file
        subject: Subject ID
        session: Session ID
        heuristic: Path to heudiconv heuristic file
        dcmconfig_json: Path to dcm2niix config JSON
        heudiconv_options: Additional heudiconv options
        output_bids_dir: Output BIDS directory
        output_info_dir: Output directory for info files
        temp_dir: Temporary directory for processing
    """
    logger.info("Single-study processing")
    logger.info("Running heudiconv on single tar file\n")

    # Find tar file
    tar_files = find_tar_files(dicoms_dir)

    if not tar_files:
        raise FileNotFoundError(f"No tar files found in {dicoms_dir}")

    if len(tar_files) > 1:
        raise ValueError(
            f"Expected single tar file but found {len(tar_files)} tar files. "
            "Use process_multi_study_heudiconv instead."
        )

    tar_file = tar_files[0]

    # Process the tar file
    info = run_heudiconv_for_study(
        tar_file,
        temp_dir,
        subject,
        session,
        heuristic,
        dcmconfig_json,
        heudiconv_options,
    )

    # Copy outputs to final locations
    output_info_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(
        info["auto_txt"], output_info_dir / f"sub-{subject}_ses-{session}_auto.txt"
    )
    shutil.copy2(
        info["dicominfo_tsv"],
        output_info_dir / f"sub-{subject}_ses-{session}_dicominfo.tsv",
    )

    # Copy BIDS directory
    output_bids_dir.mkdir(parents=True, exist_ok=True)
    if info["bids_dir"].exists():
        for src_file in info["bids_dir"].rglob("*"):
            if src_file.is_file():
                rel_path = src_file.relative_to(info["bids_dir"])
                dst_file = output_bids_dir / rel_path
                dst_file.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src_file, dst_file)

    logger.info("Single-study heudiconv processing completed successfully!")


def process_multi_study_heudiconv(
    dicoms_dir: Path,
    subject: str,
    session: str,
    heuristic: Path,
    dcmconfig_json: Path,
    heudiconv_options: str,
    output_bids_dir: Path,
    output_info_dir: Path,
    temp_dir: Path,
) -> None:
    """
    Process multiple tar files with heudiconv and merge outputs.

    Args:
        dicoms_dir: Directory containing tar files
        subject: Subject ID
        session: Session ID
        heuristic: Path to heudiconv heuristic file
        dcmconfig_json: Path to dcm2niix config JSON
        heudiconv_options: Additional heudiconv options
        output_bids_dir: Output BIDS directory
        output_info_dir: Output directory for info files
        temp_dir: Temporary directory for processing
    """
    logger.info("Multi-study processing detected")
    logger.info("Running multi-study heudiconv processing\n")

    # Find tar files
    tar_files = find_tar_files(dicoms_dir)

    if not tar_files:
        raise FileNotFoundError(f"No tar files found in {dicoms_dir}")

    if len(tar_files) == 1:
        logger.warning(
            "Only one tar file found. This function is intended for multiple studies. "
            "Consider using process_single_study_heudiconv instead."
        )

    logger.info(f"Found {len(tar_files)} tar file(s) to process")

    # Process each tar file
    info_files = []
    for tar_file in tar_files:
        info = run_heudiconv_for_study(
            tar_file,
            temp_dir,
            subject,
            session,
            heuristic,
            dcmconfig_json,
            heudiconv_options,
        )
        info_files.append(info)

    # Merge outputs
    output_info_dir.mkdir(parents=True, exist_ok=True)

    merge_auto_txt_files(
        info_files,
        output_info_dir / f"sub-{subject}_ses-{session}_auto.txt",
    )
    merge_dicominfo_files(
        info_files,
        output_info_dir / f"sub-{subject}_ses-{session}_dicominfo.tsv",
    )
    merge_bids_directories(info_files, output_bids_dir)

    logger.info("Multi-study heudiconv processing completed successfully!")
