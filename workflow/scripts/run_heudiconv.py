"""
Snakemake script to run heudiconv for DICOM to BIDS conversion.

This script handles both single-study and multi-study cases automatically.
It is called via the script: directive from the heudiconv rule.
"""

import logging
import shutil
import subprocess
import sys
from pathlib import Path

# Add workflow lib to path
sys.path.insert(0, str(Path(workflow.basedir) / "lib"))

from convert import (
    find_tar_files,
    process_multi_study_heudiconv,
)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def is_multi_study_case(dicoms_dir: Path) -> bool:
    """Check if there are multiple tar files indicating multi-study processing."""
    tar_files = find_tar_files(dicoms_dir)
    return len(tar_files) > 1


def run_single_study_heudiconv():
    """Run standard single-study heudiconv processing."""
    logger.info("Single-study processing")
    logger.info("Running standard heudiconv\n")

    cmd = (
        f"heudiconv --files {snakemake.input.dicoms_dir}"
        f" -c dcm2niix"
        f" -o {snakemake.params.out_bids}"
        f" -ss {snakemake.wildcards.session}"
        f" -s {snakemake.wildcards.subject}"
        f" -f {snakemake.input.heuristic}"
        f" --bids notop"
        f" --dcmconfig {snakemake.input.dcmconfig_json}"
        f" --overwrite"
        f" {snakemake.params.heudiconv_options}"
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error("Heudiconv failed:")
        logger.error(result.stderr)
        raise RuntimeError(f"Heudiconv failed. Check {snakemake.log[0]}")

    # Copy info files to output location
    snakemake.params.out_info_dir.mkdir(parents=True, exist_ok=True)

    in_auto_txt = (
        snakemake.params.out_bids
        / f".heudiconv/{snakemake.wildcards.subject}/ses-{snakemake.wildcards.session}/info/{snakemake.wildcards.subject}_ses-{snakemake.wildcards.session}.auto.txt"
    )
    in_dicominfo_tsv = (
        snakemake.params.out_bids
        / f".heudiconv/{snakemake.wildcards.subject}/ses-{snakemake.wildcards.session}/info/dicominfo_ses-{snakemake.wildcards.session}.tsv"
    )
    in_filegroup_json = (
        snakemake.params.out_bids
        / f".heudiconv/{snakemake.wildcards.subject}/ses-{snakemake.wildcards.session}/info/filegroup_ses-{snakemake.wildcards.session}.json"
    )

    shutil.copy2(in_auto_txt, snakemake.output.auto_txt)
    shutil.copy2(in_dicominfo_tsv, snakemake.output.dicominfo_tsv)
    shutil.copy2(in_filegroup_json, snakemake.output.filegroup_json)

    logger.info("Single-study processing completed successfully!")


def run_multi_study_heudiconv():
    """Run multi-study heudiconv processing with tar file separation and merging."""
    logger.info("Multi-study processing detected")
    logger.info("Running multi-study heudiconv processing\n")

    process_multi_study_heudiconv(
        dicoms_dir=Path(snakemake.input.dicoms_dir),
        subject=snakemake.wildcards.subject,
        session=snakemake.wildcards.session,
        heuristic=Path(snakemake.input.heuristic),
        dcmconfig_json=Path(snakemake.input.dcmconfig_json),
        heudiconv_options=snakemake.params.heudiconv_options,
        output_bids_dir=Path(snakemake.output.bids_subj_dir),
        output_info_dir=snakemake.params.out_info_dir,
        temp_dir=snakemake.params.temp_dir,
    )


def main():
    """Main entry point for the script."""
    try:
        # Check if this is a multi-study case
        if is_multi_study_case(Path(snakemake.input.dicoms_dir)):
            run_multi_study_heudiconv()
        else:
            run_single_study_heudiconv()
    except Exception as e:
        logger.error(f"Error during heudiconv processing: {e}")
        raise


if __name__ == "__main__":
    main()
