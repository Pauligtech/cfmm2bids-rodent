"""
Generate QC report for heudiconv conversion.

This script reads heudiconv metadata (*.auto.txt, dicominfo.tsv) and generates:
1. A list of series with corresponding BIDS filenames (as SVG and TSV)
2. A summary of unmapped series

Outputs are saved as SVG figures and TSV data files.
"""

import ast

import matplotlib.pyplot as plt
import pandas as pd
from lib import utils

log_file = snakemake.log[0] if snakemake.log else None
logger = utils.setup_logger(log_file)


def parse_auto_txt(auto_txt_path):
    """
    Parse the *.auto.txt file to extract BIDS mappings.

    Returns:
        dict: Mapping of series_id to BIDS path
    """

    with open(auto_txt_path) as f:
        data_str = f.read()

    data = ast.literal_eval(data_str)

    # the auto.txt file from heudiconv
    # gives mappings from bids pattern to series id
    bids_to_series = {}
    for key, series_id_list in data.items():
        bids_pattern = f"{key[0]}.{key[1][0]}"
        bids_to_series[bids_pattern] = series_id_list

    # but we want to invert this from series to bids

    # Invert mapping
    series_to_bids = {}

    for bids_str, series_list in bids_to_series.items():
        for _item, series_id in enumerate(series_list, 1):
            if isinstance(series_id, dict):
                series_id = series_id["item"]
            if series_id in series_to_bids:
                raise ValueError(
                    f"Duplicate series ID found: {series_id!r} (already mapped to {series_to_bids[series_id]!r})"
                )
            series_to_bids[series_id] = bids_str

    # Sort keys alphanumerically
    series_to_bids_sorted = dict(sorted(series_to_bids.items(), key=lambda x: x[0]))

    return series_to_bids_sorted


def load_dicominfo(dicominfo_path):
    """
    Load dicominfo.tsv file.

    Returns:
        pd.DataFrame: DICOM metadata
    """
    df = pd.read_csv(dicominfo_path, sep="\t")
    return df


def create_series_list(df, mappings, output_path, tsv_output_path=None):
    """
    Create a table showing series with BIDS filenames.

    Args:
        df: DataFrame with DICOM info
        mappings: dict mapping series_id to BIDS path
        output_path: Path to save the SVG figure
        tsv_output_path: Optional path to save the data as TSV
    """
    # Create summary DataFrame
    summary_data = []
    for _, row in df.iterrows():
        series_id = row["series_id"]
        bids_path = mappings.get(series_id, "NOT MAPPED")

        summary_data.append(
            {
                "Series ID": series_id,
                "Series Description": row["series_description"],
                "Protocol Name": row["protocol_name"],
                "Dimensions": f"{row['dim1']}×{row['dim2']}×{row['dim3']}×{row['dim4']}",
                "TR (ms)": row["TR"],
                "TE (ms)": row["TE"],
                "BIDS Path": bids_path,
            }
        )

    summary_df = pd.DataFrame(summary_data)

    # Save as TSV if path provided
    if tsv_output_path:
        summary_df.to_csv(tsv_output_path, sep="\t", index=False)
        logger.info(f"Saved series data to {tsv_output_path}")

    # Create figure with table
    fig, ax = plt.subplots(figsize=(16, max(6, len(summary_df) * 0.3)))
    ax.axis("tight")
    ax.axis("off")

    # Create table
    table = ax.table(
        cellText=summary_df.values,
        colLabels=summary_df.columns,
        cellLoc="left",
        loc="center",
        colWidths=[0.08, 0.18, 0.15, 0.12, 0.08, 0.08, 0.61],
    )

    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)

    # Style header
    for i in range(len(summary_df.columns)):
        cell = table[(0, i)]
        cell.set_facecolor("#4ECDC4")
        cell.set_text_props(weight="bold", color="white")

    # Color unmapped rows
    for i in range(len(summary_df)):
        if summary_df.iloc[i]["BIDS Path"] == "NOT MAPPED":
            for j in range(len(summary_df.columns)):
                cell = table[(i + 1, j)]
                cell.set_facecolor("#FFE5E5")

    plt.title("Series List with BIDS Mappings", fontsize=14, fontweight="bold", pad=20)
    plt.savefig(output_path, format="svg", dpi=150, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved series list to {output_path}")


def create_unmapped_summary(df, mappings, output_path):
    """
    Create a summary of unmapped series.

    Args:
        df: DataFrame with DICOM info
        mappings: dict mapping series_id to BIDS path
        output_path: Path to save the SVG figure
    """
    # Find unmapped series
    unmapped_series = []
    for _, row in df.iterrows():
        series_id = row["series_id"]
        if series_id not in mappings:
            unmapped_series.append(
                {
                    "Series ID": series_id,
                    "Series Description": row["series_description"],
                    "Protocol Name": row["protocol_name"],
                    "Files": row["series_files"],
                }
            )

    if not unmapped_series:
        # Create a simple message figure
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "All series successfully mapped to BIDS!",
            ha="center",
            va="center",
            fontsize=16,
            fontweight="bold",
            color="green",
        )
        plt.title("Unmapped Series Summary", fontsize=14, fontweight="bold", pad=20)
        plt.savefig(output_path, format="svg", dpi=150, bbox_inches="tight")
        plt.close()
        logger.info(f"Saved unmapped summary to {output_path} (all mapped)")
        return

    unmapped_df = pd.DataFrame(unmapped_series)

    # Create figure with table
    fig, ax = plt.subplots(figsize=(12, max(4, len(unmapped_df) * 0.3 + 1)))
    ax.axis("tight")
    ax.axis("off")

    # Add warning message
    warning_text = f"⚠ Warning: {len(unmapped_df)} series not mapped to BIDS"
    ax.text(
        0.5,
        0.95,
        warning_text,
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
        color="red",
        transform=ax.transAxes,
    )

    # Create table
    table = ax.table(
        cellText=unmapped_df.values,
        colLabels=unmapped_df.columns,
        cellLoc="left",
        loc="center",
        bbox=[0, 0, 1, 0.85],
    )

    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)

    # Style header
    for i in range(len(unmapped_df.columns)):
        cell = table[(0, i)]
        cell.set_facecolor("#FF6B6B")
        cell.set_text_props(weight="bold", color="white")

    plt.title("Unmapped Series Summary", fontsize=14, fontweight="bold", pad=20)
    plt.savefig(output_path, format="svg", dpi=150, bbox_inches="tight")
    plt.close()
    logger.info(f"Saved unmapped summary to {output_path}")


# main

# Construct paths to input files
dicominfo_path = snakemake.input.dicominfo_tsv
auto_txt_path = snakemake.input.auto_txt

# Load data
mappings = parse_auto_txt(auto_txt_path)
df = load_dicominfo(dicominfo_path)

series_list_output = snakemake.output.series_list
series_tsv_output = snakemake.output.series_tsv
create_series_list(df, mappings, series_list_output, series_tsv_output)

unmapped_output = snakemake.output.unmapped
create_unmapped_summary(df, mappings, unmapped_output)
