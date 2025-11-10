#!/usr/bin/env python3
"""
Script to run heudiconv on multiple tar files (multiple studies) separately,
then merge the outputs with proper series ID offsetting.

This handles the case where merge_duplicate_studies is enabled and multiple
studies are downloaded into the same directory.
"""

import argparse
import json
import logging
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


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
    study_temp_dir = temp_dir / f"study_{study_uid}"
    study_temp_dir.mkdir(parents=True, exist_ok=True)

    # Extract tar file to temp directory
    extract_dir = study_temp_dir / "dicoms"
    extract_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Extracting {tar_file.name}...")
    with tarfile.open(tar_file, "r") as tar:
        tar.extractall(extract_dir)

    # Run heudiconv on this study's DICOMs
    bids_output = study_temp_dir / "bids"
    bids_output.mkdir(parents=True, exist_ok=True)

    cmd = [
        "heudiconv",
        "--files",
        str(extract_dir),
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
    result = subprocess.run(cmd, capture_output=True, text=True)

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
    filegroup_json = heudiconv_info_dir / f"filegroup_ses-{session}.json"

    # Check that files exist
    for f in [auto_txt, dicominfo_tsv, filegroup_json]:
        if not f.exists():
            raise FileNotFoundError(f"Expected heudiconv output file not found: {f}")

    return {
        "study_uid": study_uid,
        "bids_dir": bids_output / f"sub-{subject}" / f"ses-{session}",
        "auto_txt": auto_txt,
        "dicominfo_tsv": dicominfo_tsv,
        "filegroup_json": filegroup_json,
    }


def offset_series_ids(df: pd.DataFrame, offset: int) -> pd.DataFrame:
    """Offset series_id column in dicominfo dataframe."""
    df = df.copy()
    if "series_id" in df.columns:
        df["series_id"] = df["series_id"] + offset
    return df


def merge_dicominfo_files(info_files: list[dict[str, Path]], output_tsv: Path) -> None:
    """
    Merge dicominfo.tsv files with series ID offsetting.

    Each study's series IDs are offset to avoid conflicts.
    """
    logger.info("Merging dicominfo.tsv files...")

    all_dfs = []
    offset = 0

    for i, info in enumerate(info_files):
        df = pd.read_csv(info["dicominfo_tsv"], sep="\t")

        if i > 0:
            # Calculate offset based on max series_id from previous studies
            if all_dfs:
                prev_max = max(
                    df["series_id"].max() for df in all_dfs if "series_id" in df.columns
                )
                offset = prev_max + 1000  # Add buffer to avoid conflicts
            df = offset_series_ids(df, offset)

        # Add study_uid column to track which study each series came from
        df["study_uid"] = info["study_uid"]
        all_dfs.append(df)

    merged_df = pd.concat(all_dfs, ignore_index=True)
    merged_df.to_csv(output_tsv, sep="\t", index=False)
    logger.info(f"Merged dicominfo saved to {output_tsv}")


def merge_filegroup_files(info_files: list[dict[str, Path]], output_json: Path) -> None:
    """
    Merge filegroup.json files.

    Since filegroup.json contains the heudiconv conversion info,
    we merge all studies' data together.
    """
    logger.info("Merging filegroup.json files...")

    merged_data = []

    for info in info_files:
        with open(info["filegroup_json"]) as f:
            data = json.load(f)
            # Add study_uid to each entry
            for entry in data:
                entry["study_uid"] = info["study_uid"]
            merged_data.extend(data)

    with open(output_json, "w") as f:
        json.dump(merged_data, f, indent=2)

    logger.info(f"Merged filegroup saved to {output_json}")


def merge_auto_txt_files(info_files: list[dict[str, Path]], output_txt: Path) -> None:
    """
    Merge auto.txt files.

    Simply concatenate all auto.txt files with study separators.
    """
    logger.info("Merging auto.txt files...")

    with open(output_txt, "w") as outf:
        for i, info in enumerate(info_files):
            if i > 0:
                outf.write("\n" + "=" * 80 + "\n")
                outf.write(f"Study: {info['study_uid']}\n")
                outf.write("=" * 80 + "\n\n")

            with open(info["auto_txt"]) as inf:
                outf.write(inf.read())

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

        # Copy all files from source to destination
        for src_file in src_dir.rglob("*"):
            if src_file.is_file():
                rel_path = src_file.relative_to(src_dir)
                dst_file = output_bids_dir / rel_path
                dst_file.parent.mkdir(parents=True, exist_ok=True)

                # If file already exists, add study UID to filename to avoid conflicts
                if dst_file.exists():
                    stem = dst_file.stem
                    suffix = dst_file.suffix
                    # Handle .nii.gz
                    if stem.endswith(".nii") and suffix == ".gz":
                        stem = stem[:-4]
                        suffix = ".nii.gz"
                    dst_file = (
                        dst_file.parent / f"{stem}_study-{info['study_uid']}{suffix}"
                    )

                shutil.copy2(src_file, dst_file)
                logger.debug(f"Copied {src_file} to {dst_file}")

    logger.info(f"Merged BIDS directory created at {output_bids_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Run heudiconv on multiple studies and merge outputs"
    )
    parser.add_argument(
        "--dicoms-dir",
        type=Path,
        required=True,
        help="Directory containing tar files",
    )
    parser.add_argument(
        "--subject",
        type=str,
        required=True,
        help="Subject ID",
    )
    parser.add_argument(
        "--session",
        type=str,
        required=True,
        help="Session ID",
    )
    parser.add_argument(
        "--heuristic",
        type=Path,
        required=True,
        help="Path to heudiconv heuristic file",
    )
    parser.add_argument(
        "--dcmconfig-json",
        type=Path,
        required=True,
        help="Path to dcm2niix config JSON",
    )
    parser.add_argument(
        "--heudiconv-options",
        type=str,
        default="",
        help="Additional heudiconv options",
    )
    parser.add_argument(
        "--output-bids-dir",
        type=Path,
        required=True,
        help="Output BIDS directory",
    )
    parser.add_argument(
        "--output-info-dir",
        type=Path,
        required=True,
        help="Output directory for info files",
    )
    parser.add_argument(
        "--temp-dir",
        type=Path,
        required=True,
        help="Temporary directory for processing",
    )

    args = parser.parse_args()

    # Find tar files
    tar_files = find_tar_files(args.dicoms_dir)

    if not tar_files:
        logger.error(f"No tar files found in {args.dicoms_dir}")
        sys.exit(1)

    if len(tar_files) == 1:
        logger.warning(
            "Only one tar file found. This script is intended for multiple studies. "
            "Consider using the standard heudiconv rule instead."
        )

    logger.info(f"Found {len(tar_files)} tar file(s) to process")

    # Process each tar file
    info_files = []
    for tar_file in tar_files:
        info = run_heudiconv_for_study(
            tar_file,
            args.temp_dir,
            args.subject,
            args.session,
            args.heuristic,
            args.dcmconfig_json,
            args.heudiconv_options,
        )
        info_files.append(info)

    # Merge outputs
    args.output_info_dir.mkdir(parents=True, exist_ok=True)

    merge_auto_txt_files(
        info_files,
        args.output_info_dir / f"sub-{args.subject}_ses-{args.session}_auto.txt",
    )
    merge_dicominfo_files(
        info_files,
        args.output_info_dir / f"sub-{args.subject}_ses-{args.session}_dicominfo.tsv",
    )
    merge_filegroup_files(
        info_files,
        args.output_info_dir / f"sub-{args.subject}_ses-{args.session}_filegroup.json",
    )
    merge_bids_directories(info_files, args.output_bids_dir)

    logger.info("Multi-study heudiconv processing completed successfully!")


if __name__ == "__main__":
    main()
