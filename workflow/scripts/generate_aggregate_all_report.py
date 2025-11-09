"""
Generate aggregate all HTML report with table of contents.

This script creates a top-level index/TOC page that:
1. Links to all individual subject/session reports
2. Shows summary statistics for each subject/session
3. Includes BIDS validation results for the entire dataset

The report includes:
- Overview statistics across all sessions
- Table with links to individual session reports and summary stats
- BIDS validation results (convert and fix stages)
"""

import html
import json
from pathlib import Path
from urllib.parse import quote

import pandas as pd
from lib import utils

log_file = snakemake.log[0] if snakemake.log else None
logger = utils.setup_logger(log_file)


def load_series_tsvs(series_tsv_files):
    """
    Load and aggregate all series TSV files.

    Returns:
        list: List of dicts with subject, session, and stats
    """
    logger.info(f"Loading {len(series_tsv_files)} series TSV files...")

    all_stats = []
    for tsv_file in series_tsv_files:
        subject, session = utils.extract_subject_session_from_path(tsv_file)

        if subject and session:
            df = pd.read_csv(tsv_file, sep="\t")
            n_series = len(df)
            n_unmapped = len(df[df["bids_path"] == "NOT MAPPED"])

            all_stats.append(
                {
                    "subject": subject,
                    "session": session,
                    "n_series": n_series,
                    "n_unmapped": n_unmapped,
                }
            )
        else:
            logger.warning(f"Could not extract subject/session from {tsv_file}")

    # Sort by subject then session
    all_stats.sort(key=lambda x: (x["subject"], x["session"]))

    logger.info(f"Loaded stats for {len(all_stats)} sessions")

    return all_stats


def create_html_report(session_stats, subject_reports, convert_validator, fix_validator):
    """
    Create aggregate HTML report with TOC.

    Args:
        session_stats: List of dicts with session statistics
        subject_reports: List of paths to subject report HTML files
        convert_validator: Validator JSON from convert stage
        fix_validator: Validator JSON from fix stage

    Returns:
        str: HTML content
    """
    logger.info("Creating aggregate HTML report...")
    logger.info(f"Number of session_stats: {len(session_stats)}")
    logger.info(f"Number of subject_reports: {len(subject_reports) if isinstance(subject_reports, list) else 'not a list'}")

    html_parts = []

    # Header
    html_parts.append(
        """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>cfmm2bids Aggregate Report - All Sessions</title>
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
        .report-link {
            color: #3498db;
            text-decoration: none;
            font-weight: bold;
        }
        .report-link:hover {
            text-decoration: underline;
            color: #2980b9;
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
        .warning-text {
            color: #e74c3c;
        }
        .success-text {
            color: #27ae60;
        }
    </style>
</head>
<body>
    <h1>cfmm2bids Aggregate QC Report - All Sessions</h1>
"""
    )

    # Overall Overview statistics
    html_parts.append('<div class="summary-box">')
    html_parts.append("<h2>Overall Statistics</h2>")

    if session_stats:
        n_subjects = len(set(s["subject"] for s in session_stats))
        n_sessions = len(session_stats)
        total_series = sum(s["n_series"] for s in session_stats)
        total_unmapped = sum(s["n_unmapped"] for s in session_stats)

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
        html_parts.append(f'<div class="stat-value">{total_series}</div>')
        html_parts.append("</div>")

        html_parts.append('<div class="stat">')
        html_parts.append('<div class="stat-label">Unmapped Series</div>')
        unmapped_color = "#e74c3c" if total_unmapped > 0 else "#27ae60"
        html_parts.append(
            f'<div class="stat-value" style="color: {unmapped_color};">{total_unmapped}</div>'
        )
        html_parts.append("</div>")
    else:
        html_parts.append("<p>No session data available.</p>")

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

    # Session Reports Table (TOC)
    html_parts.append('<div class="summary-box">')
    html_parts.append("<h2>Session Reports</h2>")
    html_parts.append("<p>Click on a session to view its detailed report.</p>")

    if session_stats:
        # Create mapping from (subject, session) to report filename
        # Ensure subject_reports is iterable and convert to strings if needed
        report_paths = [str(p) for p in subject_reports] if subject_reports else []
        logger.info(f"Processing {len(report_paths)} subject report paths")
        
        report_map = {}
        for report_path in report_paths:
            subject, session = utils.extract_subject_session_from_path(report_path)
            if subject and session:
                # Store the filename only
                report_map[(subject, session)] = Path(report_path).name
                logger.info(f"Added to report_map: ({subject}, {session}) -> {Path(report_path).name}")
            else:
                logger.warning(f"Could not extract subject/session from report path: {report_path}")

        logger.info(f"report_map has {len(report_map)} entries")
        if report_map:
            logger.info(f"Sample report_map keys: {list(report_map.keys())[:3]}")

        html_parts.append("<table>")
        html_parts.append("<thead><tr>")
        html_parts.append("<th>Subject</th>")
        html_parts.append("<th>Session</th>")
        html_parts.append("<th>Total Series</th>")
        html_parts.append("<th>Unmapped Series</th>")
        html_parts.append("<th>Report</th>")
        html_parts.append("</tr></thead>")

        html_parts.append("<tbody>")
        for stat in session_stats:
            subject = stat["subject"]
            session = stat["session"]
            n_series = stat["n_series"]
            n_unmapped = stat["n_unmapped"]

            # Get report filename and construct relative path with URL encoding
            report_filename = report_map.get((subject, session))
            if report_filename:
                # URL encode path components for safety
                report_link = f"sub-{quote(subject)}/ses-{quote(session)}/{quote(report_filename)}"
            else:
                # Fallback: construct expected filename even if not in map
                expected_filename = f"sub-{subject}_ses-{session}_report.html"
                report_link = f"sub-{quote(subject)}/ses-{quote(session)}/{quote(expected_filename)}"
                logger.warning(f"No report found in map for subject={subject}, session={session}, using fallback link: {report_link}")

            unmapped_class = "warning-text" if n_unmapped > 0 else "success-text"

            html_parts.append("<tr>")
            html_parts.append(f"<td>{html.escape(subject)}</td>")
            html_parts.append(f"<td>{html.escape(session)}</td>")
            html_parts.append(f"<td>{n_series}</td>")
            html_parts.append(f'<td class="{unmapped_class}">{n_unmapped}</td>')
            html_parts.append(
                f'<td><a class="report-link" href="{report_link}">View Report</a></td>'
            )
            html_parts.append("</tr>")

        html_parts.append("</tbody>")
        html_parts.append("</table>")
    else:
        html_parts.append("<p>No session reports available.</p>")

    html_parts.append("</div>")

    # Footer
    html_parts.append(
        """
</body>
</html>
"""
    )

    html_content = "\n".join(html_parts)
    logger.info("Aggregate HTML report created successfully")

    return html_content


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
logger.info("Starting aggregate all report generation...")

# Load series TSV files to get stats
series_tsv_files = snakemake.input.series_tsv
session_stats = load_series_tsvs(series_tsv_files)

# Get subject report paths - ensure it's a list
subject_reports = snakemake.input.subject_reports
if subject_reports is None:
    logger.error("subject_reports is None!")
    subject_reports = []
# If it's a single string, wrap it in a list
elif isinstance(subject_reports, str):
    logger.warning("subject_reports is a string, wrapping in list")
    subject_reports = [subject_reports]

logger.info(f"Found {len(subject_reports)} subject report files")
if len(subject_reports) > 0:
    for i, report in enumerate(subject_reports[:3]):  # Log first 3 for debugging
        logger.info(f"  Report {i+1}: {report}")
else:
    logger.error("No subject report files found!")

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
    session_stats, subject_reports, convert_validator, fix_validator
)

# Write output
output_path = snakemake.output.html_report
Path(output_path).parent.mkdir(parents=True, exist_ok=True)
with open(output_path, "w") as f:
    f.write(html_content)

logger.info(f"Aggregate all report written to {output_path}")
logger.info("Done!")
