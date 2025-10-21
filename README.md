# cfmm2bids

A Snakemake workflow for converting CFMM DICOM data to BIDS format using heudiconv.

## Features

- Download DICOM studies from CFMM
- Convert DICOM to BIDS format using heudiconv
- Generate quality control (QC) reports for each subject/session

## QC Reports

The workflow automatically generates QC reports for each subject/session after heudiconv conversion. The reports include:

1. **Gantt Chart** (`*_gantt.svg`): A timeline visualization showing all DICOM series across acquisition time, color-coded by modality (anatomical, functional, diffusion, etc.)

2. **Series List** (`*_series-list.svg`): A detailed table showing each series with:
   - Series ID and description
   - Protocol name
   - Image dimensions
   - TR and TE values
   - Corresponding BIDS filename (or "NOT MAPPED" if unmapped)

3. **Unmapped Summary** (`*_unmapped.svg`): A summary of series that were not mapped to BIDS, helping identify potential missing data or heuristic issues

QC reports are saved in the `qc/` directory with the structure:
```
qc/
└── sub-{subject}/
    └── ses-{session}/
        ├── sub-{subject}_ses-{session}_gantt.svg
        ├── sub-{subject}_ses-{session}_series-list.svg
        └── sub-{subject}_ses-{session}_unmapped.svg
```

### Manual QC Report Generation

You can also generate QC reports manually using the script:

```bash
python3 resources/generate_qc_report.py \
  --heudiconv-dir sourcedata/heudiconv/sub-001/ses-01 \
  --output-dir qc/sub-001/ses-01 \
  --subject 001 \
  --session 01
```

## Usage

1. Configure your search specifications in `config.yml`
2. Run the workflow:
   ```bash
   snakemake --cores all
   ```

## Directory Structure

```
.
├── bids/                       # BIDS-formatted output
├── sourcedata/                 # Source DICOM data
│   ├── sub-*/ses-*/           # Downloaded DICOMs
│   └── heudiconv/             # Heudiconv metadata
├── qc/                        # QC reports
├── resources/                 # Configuration files and scripts
│   ├── heuristic.py          # Heudiconv heuristic
│   ├── dcm2niix_config.json  # dcm2niix configuration
│   └── generate_qc_report.py # QC report generation script
├── config.yml                # Workflow configuration
└── Snakefile                 # Snakemake workflow
```

## Requirements

- Python 3.11+
- Snakemake
- heudiconv
- dcm2niix
- pandas
- matplotlib
- cfmm2tar

See `pixi.toml` for the complete list of dependencies.
