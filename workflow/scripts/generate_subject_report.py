"""
Generate subject-specific HTML report for a single session.

This script generates a report for one subject/session including:
1. TSV data from heudiconv QC (series.tsv file)
2. JSON provenance from post_convert_fix

The report includes:
- Overview statistics for this session
- Series table
- Summary of unmapped series
- Provenance information
"""

import html
import json
from pathlib import Path

import pandas as pd
from lib import utils

log_file = snakemake.log[0] if snakemake.log else None
logger = utils.setup_logger(log_file)


def load_series_tsv(series_tsv_file):
    """
    Load a single series TSV file.

    Returns:
        pd.DataFrame: Series data
    """
    logger.info(f"Loading series TSV file: {series_tsv_file}")

    df = pd.read_csv(series_tsv_file, sep="\t")

    logger.info(f"Loaded {len(df)} series")

    return df


def load_json_file(json_file, file_type="JSON"):
    """
    Load a single JSON file.

    Args:
        json_file: JSON file path
        file_type: Description of file type for logging

    Returns:
        dict: JSON data or None if error
    """
    logger.info(f"Loading {file_type} file: {json_file}")

    try:
        with open(json_file) as f:
            data = json.load(f)
        logger.info(f"Loaded {file_type} successfully")
        return data
    except Exception as e:
        logger.warning(f"Error loading {json_file}: {e}")
        return None


def create_html_report(subject, session, series_df, provenance_data=None):
    """
    Create HTML report for a single subject/session.

    Args:
        subject: Subject ID
        session: Session ID
        series_df: DataFrame with series data
        provenance_data: Provenance JSON data

    Returns:
        str: HTML content
    """
    logger.info(f"Creating HTML report for sub-{subject}/ses-{session}...")

    html_parts = []

    # Header
    html_parts.append(
        f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>cfmm2bids Report - sub-{html.escape(subject)} ses-{html.escape(session)}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 5px;
            margin-top: 30px;
        }}
        h3 {{
            color: #7f8c8d;
            margin-top: 20px;
        }}
        .summary-box {{
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin: 10px 0;
        }}
        .stat {{
            display: inline-block;
            margin: 10px 20px 10px 0;
            padding: 10px 15px;
            background-color: #ecf0f1;
            border-radius: 3px;
        }}
        .stat-label {{
            font-weight: bold;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
        .stat-value {{
            font-size: 1.5em;
            color: #2c3e50;
            font-weight: bold;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin: 10px 0;
        }}
        table.series-table {{
            font-size: 0.9em;
        }}
        th {{
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        td {{
            padding: 10px;
            border-bottom: 1px solid #ecf0f1;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        tr:hover {{
            background-color: #e8f4f8;
        }}
        .unmapped {{
            background-color: #ffe5e5 !important;
        }}
        pre {{
            background-color: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            font-size: 0.9em;
        }}
        .json-viewer {{
            background-color: white;
            padding: 10px;
            border-radius: 5px;
            border: 1px solid #ddd;
            margin: 10px 0;
            max-height: 400px;
            overflow-y: auto;
        }}
        details {{
            background-color: white;
            margin: 10px 0;
            border-radius: 3px;
            overflow: hidden;
        }}
        summary {{
            background-color: #3498db;
            color: white;
            cursor: pointer;
            padding: 10px;
            font-size: 1em;
            font-weight: bold;
            list-style: none;
            user-select: none;
        }}
        summary:hover {{
            background-color: #2980b9;
        }}
        summary::-webkit-details-marker {{
            display: none;
        }}
        summary::before {{
            content: '▶';
            display: inline-block;
            margin-right: 8px;
            transition: transform 0.2s;
        }}
        details[open] > summary::before {{
            transform: rotate(90deg);
        }}
        details > div {{
            padding: 10px;
        }}
    </style>
</head>
<body>
    <h1>cfmm2bids QC Report</h1>
    <h2>Subject: {html.escape(subject)} | Session: {html.escape(session)}</h2>
"""
    )

    # Overview statistics
    html_parts.append('<div class="summary-box">')
    html_parts.append("<h2>Overview Statistics</h2>")

    if not series_df.empty:
        n_series = len(series_df)
        n_unmapped = len(series_df[series_df["bids_path"] == "NOT MAPPED"])

        html_parts.append('<div class="stat">')
        html_parts.append('<div class="stat-label">Total Series</div>')
        html_parts.append(f'<div class="stat-value">{n_series}</div>')
        html_parts.append("</div>")

        html_parts.append('<div class="stat">')
        html_parts.append('<div class="stat-label">Unmapped Series</div>')
        unmapped_color = "#e74c3c" if n_unmapped > 0 else "#27ae60"
        html_parts.append(
            f'<div class="stat-value" style="color: {unmapped_color};">{n_unmapped}</div>'
        )
        html_parts.append("</div>")
    else:
        html_parts.append("<p>No series data available.</p>")

    html_parts.append("</div>")

    # Series Table
    html_parts.append('<div class="summary-box">')
    html_parts.append("<h2>Series Table</h2>")

    if not series_df.empty:
        html_parts.append('<div style="overflow-x: auto;">')
        html_parts.append(series_df_to_html_table(series_df))
        html_parts.append("</div>")
    else:
        html_parts.append("<p>No series data available.</p>")

    html_parts.append("</div>")

    # Provenance Information
    if provenance_data:
        html_parts.append('<div class="summary-box">')
        html_parts.append("<h2>Post-Conversion Fix Provenance</h2>")

        html_parts.append("<details>")
        html_parts.append("<summary>View Provenance Data</summary>")
        html_parts.append("<div>")
        html_parts.append('<div class="json-viewer">')
        html_parts.append(
            f"<pre>{html.escape(json.dumps(provenance_data, indent=2))}</pre>"
        )
        html_parts.append("</div>")
        html_parts.append("</div>")
        html_parts.append("</details>")

        html_parts.append("</div>")

    # Footer
    html_parts.append(
        """
</body>
</html>
"""
    )

    html_content = "\n".join(html_parts)
    logger.info("HTML report created successfully")

    return html_content


def series_df_to_html_table(df):
    """
    Convert series DataFrame to HTML table with custom styling.

    Args:
        df: DataFrame with series data

    Returns:
        str: HTML table
    """
    html_parts = ['<table class="series-table">']

    # Header
    html_parts.append("<thead><tr>")
    for col in df.columns:
        html_parts.append(f"<th>{html.escape(col)}</th>")
    html_parts.append("</tr></thead>")

    # Body
    html_parts.append("<tbody>")
    for _, row in df.iterrows():
        # Check if unmapped
        is_unmapped = row.get("bids_path") == "NOT MAPPED"
        row_class = ' class="unmapped"' if is_unmapped else ""

        html_parts.append(f"<tr{row_class}>")
        for col in df.columns:
            value = row[col]
            # Format value
            value = "" if pd.isna(value) else str(value)
            # Escape HTML to prevent XSS
            value = html.escape(value)

            html_parts.append(f"<td>{value}</td>")
        html_parts.append("</tr>")
    html_parts.append("</tbody>")

    html_parts.append("</table>")

    return "\n".join(html_parts)


# Main execution
logger.info("Starting subject report generation...")

# Extract subject and session from wildcards
subject = snakemake.wildcards.subject
session = snakemake.wildcards.session

logger.info(f"Generating report for sub-{subject}/ses-{session}")

# Load series TSV file
series_tsv_file = snakemake.input.series_tsv
series_df = load_series_tsv(series_tsv_file)

# Load provenance JSON file
prov_json = snakemake.input.get("provenance_json", False)
if prov_json:
    provenance_data = load_json_file(prov_json, "provenance JSON")

    # Create HTML report
    html_content = create_html_report(subject, session, series_df, provenance_data)
else:
    # Create HTML report
    html_content = create_html_report(subject, session, series_df)


# Write output
output_path = snakemake.output.html_report
Path(output_path).parent.mkdir(parents=True, exist_ok=True)
with open(output_path, "w") as f:
    f.write(html_content)

logger.info(f"Subject report written to {output_path}")
logger.info("Done!")
