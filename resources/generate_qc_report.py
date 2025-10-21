#!/usr/bin/env python3
"""
Generate QC report for heudiconv conversion.

This script reads heudiconv metadata (*.auto.txt, dicominfo.tsv) and generates:
1. A Gantt-style chart showing all series across acquisition time
2. A list of series with corresponding BIDS filenames
3. A summary of unmapped series

Outputs are saved as SVG figures.
"""

import argparse
import re
import sys
from pathlib import Path
from datetime import datetime, timedelta

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.dates import DateFormatter, AutoDateLocator


def parse_auto_txt(auto_txt_path):
    """
    Parse the *.auto.txt file to extract BIDS mappings.
    
    Returns:
        dict: Mapping of series_id to BIDS path
    """
    mappings = {}
    
    with open(auto_txt_path, 'r') as f:
        lines = f.readlines()
    
    current_series_id = None
    for line in lines:
        line = line.strip()
        
        # Match series info line: "seqinfo: <series_id> <description> [...]"
        series_match = re.match(r'seqinfo:\s+(\d+)\s+(.+?)(?:\s+\[.*\])?$', line)
        if series_match:
            current_series_id = int(series_match.group(1))
            continue
        
        # Match BIDS mapping line: "  > <bids_path>"
        bids_match = re.match(r'\s*>\s+(.+)$', line)
        if bids_match and current_series_id is not None:
            bids_path = bids_match.group(1)
            mappings[current_series_id] = bids_path
            current_series_id = None
    
    return mappings


def load_dicominfo(dicominfo_path):
    """
    Load dicominfo.tsv file.
    
    Returns:
        pd.DataFrame: DICOM metadata
    """
    df = pd.read_csv(dicominfo_path, sep='\t')
    return df


def calculate_acquisition_times(df):
    """
    Calculate start and end acquisition times for each series.
    
    Args:
        df: DataFrame with DICOM info
        
    Returns:
        pd.DataFrame: DataFrame with start_time and end_time columns
    """
    # Parse date and time
    df['datetime'] = pd.to_datetime(
        df['date'].astype(str) + df['time'].astype(str).str.zfill(6),
        format='%Y%m%d%H%M%S',
        errors='coerce'
    )
    
    # Calculate duration based on TR and number of volumes (dim4)
    # For 3D images (dim4=1), duration is essentially instantaneous (use TR as proxy)
    # For 4D images, duration = TR * dim4 / 1000 (TR is in ms)
    df['duration_seconds'] = (df['TR'] * df['dim4']) / 1000.0
    
    # Calculate end time
    df['end_time'] = df['datetime'] + pd.to_timedelta(df['duration_seconds'], unit='s')
    
    return df


def create_gantt_chart(df, mappings, output_path):
    """
    Create a Gantt-style chart showing series across time.
    
    Args:
        df: DataFrame with DICOM info including datetime and end_time
        mappings: dict mapping series_id to BIDS path
        output_path: Path to save the SVG figure
    """
    fig, ax = plt.subplots(figsize=(14, max(8, len(df) * 0.4)))
    
    # Color mapping for different modalities
    modality_colors = {
        'anat': '#FF6B6B',
        'func': '#4ECDC4',
        'dwi': '#95E1D3',
        'fmap': '#F38181',
        'other': '#CCCCCC',
        'unmapped': '#999999',
    }
    
    y_positions = []
    y_labels = []
    
    for idx, row in df.iterrows():
        series_id = row['series_id']
        series_desc = row['series_description']
        start_time = row['datetime']
        end_time = row['end_time']
        
        # Check if series is mapped
        if series_id in mappings:
            bids_path = mappings[series_id]
            # Determine modality from BIDS path
            if '/anat/' in bids_path:
                color = modality_colors['anat']
                modality = 'anat'
            elif '/func/' in bids_path:
                color = modality_colors['func']
                modality = 'func'
            elif '/dwi/' in bids_path:
                color = modality_colors['dwi']
                modality = 'dwi'
            elif '/fmap/' in bids_path:
                color = modality_colors['fmap']
                modality = 'fmap'
            else:
                color = modality_colors['other']
                modality = 'other'
            label = f"[{series_id}] {series_desc}"
        else:
            color = modality_colors['unmapped']
            modality = 'unmapped'
            label = f"[{series_id}] {series_desc} (unmapped)"
        
        # Skip if datetime is NaT
        if pd.isna(start_time) or pd.isna(end_time):
            continue
        
        # Create bar
        duration = (end_time - start_time).total_seconds() / 60.0  # Convert to minutes
        ax.barh(idx, duration, left=start_time, height=0.7, 
                color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        y_positions.append(idx)
        y_labels.append(label)
    
    # Format axes
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel('Acquisition Time', fontsize=12)
    ax.set_ylabel('Series', fontsize=12)
    ax.set_title('DICOM Acquisition Timeline (Gantt Chart)', fontsize=14, fontweight='bold')
    
    # Format x-axis to show time
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    # Use AutoDateLocator instead of MinuteLocator to avoid too many ticks
    from matplotlib.dates import AutoDateLocator
    ax.xaxis.set_major_locator(AutoDateLocator())
    
    # Add legend
    legend_elements = [
        mpatches.Patch(color=modality_colors['anat'], label='Anatomical'),
        mpatches.Patch(color=modality_colors['func'], label='Functional'),
        mpatches.Patch(color=modality_colors['dwi'], label='Diffusion'),
        mpatches.Patch(color=modality_colors['fmap'], label='Fieldmap'),
        mpatches.Patch(color=modality_colors['other'], label='Other'),
        mpatches.Patch(color=modality_colors['unmapped'], label='Unmapped'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    # Adjust layout
    plt.tight_layout()
    plt.grid(axis='x', alpha=0.3)
    
    # Save figure
    plt.savefig(output_path, format='svg', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved Gantt chart to {output_path}")


def create_series_list(df, mappings, output_path):
    """
    Create a table showing series with BIDS filenames.
    
    Args:
        df: DataFrame with DICOM info
        mappings: dict mapping series_id to BIDS path
        output_path: Path to save the SVG figure
    """
    # Create summary DataFrame
    summary_data = []
    for _, row in df.iterrows():
        series_id = row['series_id']
        bids_path = mappings.get(series_id, 'NOT MAPPED')
        
        summary_data.append({
            'Series ID': series_id,
            'Series Description': row['series_description'],
            'Protocol Name': row['protocol_name'],
            'Dimensions': f"{row['dim1']}×{row['dim2']}×{row['dim3']}×{row['dim4']}",
            'TR (ms)': row['TR'],
            'TE (ms)': row['TE'],
            'BIDS Path': bids_path,
        })
    
    summary_df = pd.DataFrame(summary_data)
    
    # Create figure with table
    fig, ax = plt.subplots(figsize=(16, max(6, len(summary_df) * 0.3)))
    ax.axis('tight')
    ax.axis('off')
    
    # Create table
    table = ax.table(cellText=summary_df.values,
                     colLabels=summary_df.columns,
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.08, 0.18, 0.15, 0.12, 0.08, 0.08, 0.31])
    
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)
    
    # Style header
    for i in range(len(summary_df.columns)):
        cell = table[(0, i)]
        cell.set_facecolor('#4ECDC4')
        cell.set_text_props(weight='bold', color='white')
    
    # Color unmapped rows
    for i in range(len(summary_df)):
        if summary_df.iloc[i]['BIDS Path'] == 'NOT MAPPED':
            for j in range(len(summary_df.columns)):
                cell = table[(i + 1, j)]
                cell.set_facecolor('#FFE5E5')
    
    plt.title('Series List with BIDS Mappings', fontsize=14, fontweight='bold', pad=20)
    plt.savefig(output_path, format='svg', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved series list to {output_path}")


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
        series_id = row['series_id']
        if series_id not in mappings:
            unmapped_series.append({
                'Series ID': series_id,
                'Series Description': row['series_description'],
                'Protocol Name': row['protocol_name'],
                'Files': row['series_files'],
            })
    
    if not unmapped_series:
        # Create a simple message figure
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.axis('off')
        ax.text(0.5, 0.5, 'All series successfully mapped to BIDS!',
                ha='center', va='center', fontsize=16, fontweight='bold',
                color='green')
        plt.title('Unmapped Series Summary', fontsize=14, fontweight='bold', pad=20)
        plt.savefig(output_path, format='svg', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved unmapped summary to {output_path} (all mapped)")
        return
    
    unmapped_df = pd.DataFrame(unmapped_series)
    
    # Create figure with table
    fig, ax = plt.subplots(figsize=(12, max(4, len(unmapped_df) * 0.3 + 1)))
    ax.axis('tight')
    ax.axis('off')
    
    # Add warning message
    warning_text = f"⚠ Warning: {len(unmapped_df)} series not mapped to BIDS"
    ax.text(0.5, 0.95, warning_text,
            ha='center', va='top', fontsize=12, fontweight='bold',
            color='red', transform=ax.transAxes)
    
    # Create table
    table = ax.table(cellText=unmapped_df.values,
                     colLabels=unmapped_df.columns,
                     cellLoc='left',
                     loc='center',
                     bbox=[0, 0, 1, 0.85])
    
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.8)
    
    # Style header
    for i in range(len(unmapped_df.columns)):
        cell = table[(0, i)]
        cell.set_facecolor('#FF6B6B')
        cell.set_text_props(weight='bold', color='white')
    
    plt.title('Unmapped Series Summary', fontsize=14, fontweight='bold', pad=20)
    plt.savefig(output_path, format='svg', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved unmapped summary to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate QC report for heudiconv conversion'
    )
    parser.add_argument(
        '--heudiconv-dir',
        required=True,
        type=Path,
        help='Path to heudiconv directory (e.g., sourcedata/heudiconv/sub-001/ses-01)'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        type=Path,
        help='Directory to save QC reports'
    )
    parser.add_argument(
        '--subject',
        required=True,
        help='Subject ID'
    )
    parser.add_argument(
        '--session',
        required=True,
        help='Session ID'
    )
    
    args = parser.parse_args()
    
    # Construct paths to input files
    info_dir = args.heudiconv_dir / 'info'
    dicominfo_path = info_dir / 'dicominfo.tsv'
    auto_txt_pattern = info_dir / f'sub-{args.subject}_ses-{args.session}.auto.txt'
    
    # Find auto.txt file (may have different naming)
    auto_txt_files = list(info_dir.glob('*.auto.txt'))
    if not auto_txt_files:
        print(f"Error: No *.auto.txt file found in {info_dir}")
        sys.exit(1)
    auto_txt_path = auto_txt_files[0]
    
    # Check if files exist
    if not dicominfo_path.exists():
        print(f"Error: dicominfo.tsv not found at {dicominfo_path}")
        sys.exit(1)
    
    print(f"Reading heudiconv metadata from {args.heudiconv_dir}")
    print(f"  - dicominfo.tsv: {dicominfo_path}")
    print(f"  - auto.txt: {auto_txt_path}")
    
    # Load data
    mappings = parse_auto_txt(auto_txt_path)
    df = load_dicominfo(dicominfo_path)
    df = calculate_acquisition_times(df)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate reports
    print(f"\nGenerating QC reports for sub-{args.subject} ses-{args.session}...")
    
    gantt_output = args.output_dir / f'sub-{args.subject}_ses-{args.session}_gantt.svg'
    create_gantt_chart(df, mappings, gantt_output)
    
    series_list_output = args.output_dir / f'sub-{args.subject}_ses-{args.session}_series-list.svg'
    create_series_list(df, mappings, series_list_output)
    
    unmapped_output = args.output_dir / f'sub-{args.subject}_ses-{args.session}_unmapped.svg'
    create_unmapped_summary(df, mappings, unmapped_output)
    
    print(f"\n✓ QC reports generated successfully in {args.output_dir}")
    print(f"  - Gantt chart: {gantt_output.name}")
    print(f"  - Series list: {series_list_output.name}")
    print(f"  - Unmapped summary: {unmapped_output.name}")


if __name__ == '__main__':
    main()
