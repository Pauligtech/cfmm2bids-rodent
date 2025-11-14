"""
Snakemake script to run heudiconv for DICOM to BIDS conversion.

This script handles both single-study and multi-study cases automatically.
It is called via the script: directive from the heudiconv rule.
"""

from pathlib import Path

from lib.convert import (
    find_tar_files,
    process_multi_study_heudiconv,
    process_single_study_heudiconv,
)
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)


def is_multi_study_case(dicoms_dir: Path) -> bool:
    """Check if there are multiple tar files indicating multi-study processing."""
    tar_files = find_tar_files(dicoms_dir)
    return len(tar_files) > 1


def main():
    """Main entry point for the script."""
    try:
        # Create bids-temp directory for processing
        temp_dir = Path("bids-temp")
        temp_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using temp directory: {temp_dir.absolute()}")

        # Check if this is a multi-study case
        if is_multi_study_case(Path(snakemake.input.dicoms_dir)):
            process_multi_study_heudiconv(
                dicoms_dir=Path(snakemake.input.dicoms_dir),
                subject=snakemake.wildcards.subject,
                session=snakemake.wildcards.session,
                heuristic=Path(snakemake.input.heuristic),
                dcmconfig_json=Path(snakemake.input.dcmconfig_json),
                heudiconv_options=snakemake.params.heudiconv_options,
                output_bids_dir=Path(snakemake.output.bids_subj_dir),
                output_info_dir=snakemake.params.out_info_dir,
                temp_dir=temp_dir,
            )
        else:
            process_single_study_heudiconv(
                dicoms_dir=Path(snakemake.input.dicoms_dir),
                subject=snakemake.wildcards.subject,
                session=snakemake.wildcards.session,
                heuristic=Path(snakemake.input.heuristic),
                dcmconfig_json=Path(snakemake.input.dcmconfig_json),
                heudiconv_options=snakemake.params.heudiconv_options,
                output_bids_dir=Path(snakemake.output.bids_subj_dir),
                output_info_dir=snakemake.params.out_info_dir,
                temp_dir=temp_dir,
            )
    except Exception as e:
        logger.error(f"Error during heudiconv processing: {e}")
        raise


if __name__ == "__main__":
    main()
