"""
Generate aggregate HTML report for all sessions.

This script aggregates:
1. TSV data from heudiconv QC (series.tsv files)
2. JSON metadata from heudiconv (filegroup.json)
3. JSON provenance from post_convert_fix
4. Validator results (JSON from convert and fix stages)

The report includes:
- Overview statistics
- Aggregated series table (sorted by subject, session)
- Summary of unmapped series
- Validation results
- Provenance information
"""

import html
import json
from pathlib import Path

import pandas as pd
from lib import utils

log_file = snakemake.log[0] if snakemake.log else None
logger = utils.setup_logger(log_file)


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


def load_series_tsvs(series_tsv_files):
    """
    Load and aggregate all series TSV files.

    Returns:
        pd.DataFrame: Aggregated series data with subject and session columns
    """
    logger.info(f"Loading {len(series_tsv_files)} series TSV files...")

    all_data = []
    for tsv_file in series_tsv_files:
        subject, session = extract_subject_session_from_path(tsv_file)

        if subject and session:
            df = pd.read_csv(tsv_file, sep="\t")
            df.insert(0, "subject", subject)
            df.insert(1, "session", session)
            all_data.append(df)
        else:
            logger.warning(f"Could not extract subject/session from {tsv_file}")

    if not all_data:
        logger.warning("No series TSV files loaded")
        return pd.DataFrame()

    # Concatenate all dataframes
    combined_df = pd.concat(all_data, ignore_index=True)

    # Sort by subject then session
    combined_df = combined_df.sort_values(by=["subject", "session"]).reset_index(
        drop=True
    )

    logger.info(
        f"Loaded {len(combined_df)} series across {len(series_tsv_files)} sessions"
    )

    return combined_df


def load_json_files(json_files, file_type="JSON"):
    """
    Load multiple JSON files into a list of dicts.

    Args:
        json_files: List of JSON file paths
        file_type: Description of file type for logging

    Returns:
        list: List of dicts, each containing the JSON data plus metadata
    """
    logger.info(f"Loading {len(json_files)} {file_type} files...")

    all_data = []
    for json_file in json_files:
        subject, session = extract_subject_session_from_path(json_file)

        try:
            with open(json_file) as f:
                data = json.load(f)

            all_data.append(
                {
                    "subject": subject,
                    "session": session,
                    "file": str(json_file),
                    "data": data,
                }
            )
        except Exception as e:
            logger.warning(f"Error loading {json_file}: {e}")

    # Sort by subject then session
    all_data.sort(key=lambda x: (x.get("subject", ""), x.get("session", "")))

    logger.info(f"Loaded {len(all_data)} {file_type} files")

    return all_data


def create_html_report(
    series_df, filegroup_data, provenance_data, convert_validator, fix_validator
):
    """
    Create HTML report from aggregated data.

    Args:
        series_df: DataFrame with aggregated series data
        filegroup_data: List of filegroup JSON data
        provenance_data: List of provenance JSON data
        convert_validator: Validator JSON from convert stage
        fix_validator: Validator JSON from fix stage

    Returns:
        str: HTML content
    """
    logger.info("Creating HTML report...")

    html_parts = []

    # Header
    html_parts.append(
        """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>cfmm2bids Aggregate Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 5px;
            margin-top: 30px;
        }
        h3 {
            color: #7f8c8d;
            margin-top: 20px;
        }
        .summary-box {
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin: 10px 0;
        }
        .stat {
            display: inline-block;
            margin: 10px 20px 10px 0;
            padding: 10px 15px;
            background-color: #ecf0f1;
            border-radius: 3px;
        }
        .stat-label {
            font-weight: bold;
            color: #7f8c8d;
            font-size: 0.9em;
        }
        .stat-value {
            font-size: 1.5em;
            color: #2c3e50;
            font-weight: bold;
        }
        table {
            border-collapse: collapse;
            width: 100%;
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin: 10px 0;
        }
        table.series-table {
            font-size: 0.9em;
        }
        th {
            background-color: #3498db;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }
        td {
            padding: 10px;
            border-bottom: 1px solid #ecf0f1;
        }
        tr:nth-child(even) {
            background-color: #f8f9fa;
        }
        tr:hover {
            background-color: #e8f4f8;
        }
        .unmapped {
            background-color: #ffe5e5 !important;
        }
        .validation-pass {
            color: #27ae60;
            font-weight: bold;
        }
        .validation-fail {
            color: #e74c3c;
            font-weight: bold;
        }
        .validation-warn {
            color: #f39c12;
            font-weight: bold;
        }
        pre {
            background-color: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            font-size: 0.9em;
        }
        .json-viewer {
            background-color: white;
            padding: 10px;
            border-radius: 5px;
            border: 1px solid #ddd;
            margin: 10px 0;
            max-height: 400px;
            overflow-y: auto;
        }
        .collapsible {
            background-color: #3498db;
            color: white;
            cursor: pointer;
            padding: 10px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            font-size: 1em;
            font-weight: bold;
            margin-top: 10px;
            border-radius: 3px;
        }
        .collapsible:hover {
            background-color: #2980b9;
        }
        .collapsible:after {
            content: '\\002B'; /* Plus sign */
            font-weight: bold;
            float: right;
        }
        .collapsible.active:after {
            content: '\\2212'; /* Minus sign */
        }
        .content {
            padding: 0 10px;
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.2s ease-out;
            background-color: white;
        }
    </style>
</head>
<body>
    <h1>cfmm2bids Aggregate QC Report</h1>
"""
    )

    # Overview statistics
    html_parts.append('<div class="summary-box">')
    html_parts.append("<h2>Overview Statistics</h2>")

    if not series_df.empty:
        n_subjects = series_df["subject"].nunique()
        n_sessions = len(series_df.groupby(["subject", "session"]))
        n_series = len(series_df)
        n_unmapped = len(series_df[series_df["bids_path"] == "NOT MAPPED"])

        html_parts.append('<div class="stat">')
        html_parts.append('<div class="stat-label">Subjects</div>')
        html_parts.append(f'<div class="stat-value">{n_subjects}</div>')
        html_parts.append("</div>")

        html_parts.append('<div class="stat">')
        html_parts.append('<div class="stat-label">Sessions</div>')
        html_parts.append(f'<div class="stat-value">{n_sessions}</div>')
        html_parts.append("</div>")

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

    # BIDS Validation Results
    html_parts.append('<div class="summary-box">')
    html_parts.append("<h2>BIDS Validation Results</h2>")

    # Convert stage validation
    html_parts.append("<h3>Convert Stage Validation</h3>")
    if convert_validator:
        html_parts.append(format_validator_summary(convert_validator))
    else:
        html_parts.append("<p>No validation data available for convert stage.</p>")

    # Fix stage validation
    html_parts.append("<h3>Fix Stage Validation</h3>")
    if fix_validator:
        html_parts.append(format_validator_summary(fix_validator))
    else:
        html_parts.append("<p>No validation data available for fix stage.</p>")

    html_parts.append("</div>")

    # Aggregated Series Table
    html_parts.append('<div class="summary-box">')
    html_parts.append("<h2>All Series (Sorted by Subject, Session)</h2>")

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

        for prov in provenance_data:
            subj = prov.get("subject", "unknown")
            sess = prov.get("session", "unknown")
            data = prov.get("data", {})

            html_parts.append(
                f'<button class="collapsible">sub-{html.escape(subj)} / ses-{html.escape(sess)}</button>'
            )
            html_parts.append('<div class="content">')
            html_parts.append('<div class="json-viewer">')
            html_parts.append(f"<pre>{html.escape(json.dumps(data, indent=2))}</pre>")
            html_parts.append("</div>")
            html_parts.append("</div>")

        html_parts.append("</div>")

    # Filegroup Information (collapsible)
    if filegroup_data:
        html_parts.append('<div class="summary-box">')
        html_parts.append("<h2>Heudiconv Filegroup Metadata</h2>")
        html_parts.append(
            "<p>Click to expand filegroup data for each subject/session:</p>"
        )

        for fg in filegroup_data:
            subj = fg.get("subject", "unknown")
            sess = fg.get("session", "unknown")
            data = fg.get("data", {})

            html_parts.append(
                f'<button class="collapsible">sub-{html.escape(subj)} / ses-{html.escape(sess)}</button>'
            )
            html_parts.append('<div class="content">')
            html_parts.append('<div class="json-viewer">')
            html_parts.append(f"<pre>{html.escape(json.dumps(data, indent=2))}</pre>")
            html_parts.append("</div>")
            html_parts.append("</div>")

        html_parts.append("</div>")

    # Footer with JavaScript for collapsibles
    html_parts.append(
        """
    <script>
    var coll = document.getElementsByClassName("collapsible");
    var i;
    for (i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.maxHeight){
                content.style.maxHeight = null;
            } else {
                content.style.maxHeight = content.scrollHeight + "px";
            }
        });
    }
    </script>
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


def format_validator_summary(validator_data):
    """
    Format validator JSON data into HTML summary.

    Supports the bids-validator-deno JSON format with issues.issues array.

    Args:
        validator_data: Dictionary with validator results

    Returns:
        str: HTML summary
    """
    from collections import defaultdict

    html_parts = []

    # Parse the validator data structure
    issues_data = validator_data.get("issues", {})
    issues_list = issues_data.get("issues", [])
    code_messages = issues_data.get("codeMessages", {})
    summary = validator_data.get("summary", {})

    # Separate errors and warnings by severity
    errors = [issue for issue in issues_list if issue.get("severity") == "error"]
    warnings = [issue for issue in issues_list if issue.get("severity") == "warning"]

    # Check if valid (no errors means valid)
    is_valid = len(errors) == 0

    if is_valid:
        html_parts.append('<p class="validation-pass">✓ Dataset is BIDS compliant</p>')
    else:
        html_parts.append(
            '<p class="validation-fail">✗ Dataset has validation issues</p>'
        )

    # Display summary section
    html_parts.append("<h4>Summary:</h4>")
    html_parts.append('<div class="json-viewer" style="max-height: 300px;">')
    html_parts.append(f"<pre>{html.escape(json.dumps(summary, indent=2))}</pre>")
    html_parts.append("</div>")

    # Display counts
    html_parts.append(
        f'<p><strong>Errors:</strong> <span class="validation-fail">{len(errors)}</span></p>'
    )
    html_parts.append(
        f'<p><strong>Warnings:</strong> <span class="validation-warn">{len(warnings)}</span></p>'
    )

    # Helper function to format issues hierarchically
    def format_issues_hierarchy(issues, severity_label):
        if not issues:
            return

        html_parts.append(f"<h4>{severity_label}:</h4>")

        # Organize issues: code -> subCode -> locations
        hierarchy = defaultdict(lambda: defaultdict(list))

        for issue in issues:
            code = issue.get("code", "Unknown")
            sub_code = issue.get("subCode", None)
            location = issue.get("location", "Unknown")
            issue_msg = issue.get("issueMessage", "")

            hierarchy[code][sub_code].append(
                {"location": location, "issueMessage": issue_msg}
            )

        # Generate collapsible HTML for each code
        for code in sorted(hierarchy.keys()):
            code_msg = code_messages.get(code, "")
            subcodes = hierarchy[code]

            # Create a unique ID for this collapsible section
            section_id = f"{severity_label.lower()}_{code}_{id(code)}"

            html_parts.append(
                f'<button class="collapsible" data-target="{section_id}">'
            )
            html_parts.append(f"{html.escape(code)} ({sum(len(locs) for locs in subcodes.values())} issues)")
            html_parts.append("</button>")
            html_parts.append(f'<div class="content" id="{section_id}">')

            # Show code-level message if available
            if code_msg:
                html_parts.append(f"<p><em>{html.escape(code_msg)}</em></p>")

            # Show subcodes (or locations if no subcode)
            for sub_code in sorted(subcodes.keys()):
                locations = subcodes[sub_code]

                if sub_code is not None:
                    # Has subCode - create another level of hierarchy
                    subcode_id = f"{section_id}_{sub_code}_{id(sub_code)}"
                    html_parts.append(
                        f'<button class="collapsible" style="margin-left: 20px;" data-target="{subcode_id}">'
                    )
                    html_parts.append(
                        f"SubCode: {html.escape(sub_code)} ({len(locations)} locations)"
                    )
                    html_parts.append("</button>")
                    html_parts.append(f'<div class="content" id="{subcode_id}">')

                # Show locations
                html_parts.append('<ul style="margin-left: 20px;">')
                for loc_info in locations:
                    location = loc_info["location"]
                    issue_msg = loc_info["issueMessage"]

                    html_parts.append(f"<li><strong>{html.escape(location)}</strong>")
                    if issue_msg:
                        html_parts.append(
                            f"<br/><em>{html.escape(issue_msg.strip())}</em>"
                        )
                    html_parts.append("</li>")
                html_parts.append("</ul>")

                if sub_code is not None:
                    html_parts.append("</div>")  # Close subcode content

            html_parts.append("</div>")  # Close code content

    # Format errors
    format_issues_hierarchy(errors, "Errors")

    # Format warnings
    format_issues_hierarchy(warnings, "Warnings")

    return "\n".join(html_parts)


# Main execution
logger.info("Starting aggregate report generation...")

# Load series TSV files
series_tsv_files = snakemake.input.series_tsv
series_df = load_series_tsvs(series_tsv_files)

# Load filegroup JSON files
filegroup_files = snakemake.input.filegroup_json
filegroup_data = load_json_files(filegroup_files, "filegroup JSON")

# Load provenance JSON files
provenance_files = snakemake.input.provenance_json
provenance_data = load_json_files(provenance_files, "provenance JSON")

# Load validator JSON files
convert_validator = None
fix_validator = None

try:
    with open(snakemake.input.convert_validator_json) as f:
        convert_validator = json.load(f)
    logger.info("Loaded convert stage validator JSON")
except Exception as e:
    logger.warning(f"Could not load convert validator JSON: {e}")

try:
    with open(snakemake.input.fix_validator_json) as f:
        fix_validator = json.load(f)
    logger.info("Loaded fix stage validator JSON")
except Exception as e:
    logger.warning(f"Could not load fix validator JSON: {e}")

# Create HTML report
html_content = create_html_report(
    series_df, filegroup_data, provenance_data, convert_validator, fix_validator
)

# Write output
output_path = snakemake.output.html_report
Path(output_path).parent.mkdir(parents=True, exist_ok=True)
with open(output_path, "w") as f:
    f.write(html_content)

logger.info(f"Aggregate report written to {output_path}")
logger.info("Done!")
