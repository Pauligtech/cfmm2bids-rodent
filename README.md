# cfmm2bids

A Snakemake workflow for converting CFMM DICOM data to BIDS format using heudiconv.

## Features

- Download DICOM studies from CFMM
- Convert DICOM to BIDS format using heudiconv
- Generate quality control (QC) reports for each subject/session

## QC Reports

The workflow automatically generates QC reports for each subject/session after heudiconv conversion. The reports include:

1. **Series List** (`*_series.svg`): A detailed table showing each series with:
   - Series ID and description
   - Protocol name
   - Image dimensions
   - TR and TE values
   - Corresponding BIDS filename (or "NOT MAPPED" if unmapped)

2. **Unmapped Summary** (`*_unmapped.svg`): A summary of series that were not mapped to BIDS, helping identify potential missing data or heuristic issues

QC reports are saved in the `sourcedata/qc/` directory with the structure:
```
sourcedata/qc/
└── sub-{subject}/
    └── ses-{session}/
        ├── sub-{subject}_ses-{session}_series.svg
        └── sub-{subject}_ses-{session}_unmapped.svg
```

**Note:** The QC report generation is integrated into the Snakemake workflow as a script directive and cannot be run manually as a standalone CLI tool.

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
│   ├── heudiconv/             # Heudiconv metadata
│   └── qc/                    # QC reports
├── workflow/                  # Workflow files
│   ├── Snakefile              # Snakemake workflow
│   ├── lib/                   # Python modules
│   │   └── query_filter.py   # DICOM query filtering
│   └── scripts/               # Workflow scripts
│       └── generate_qc_report.py # QC report generation script
├── resources/                 # Configuration files
│   ├── heuristic.py          # Heudiconv heuristic
│   └── dcm2niix_config.json  # dcm2niix configuration
└── config.yml                # Workflow configuration
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
