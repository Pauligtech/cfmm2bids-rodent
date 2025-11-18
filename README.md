# cfmm2bids

A Snakemake workflow for converting CFMM DICOM data to BIDS format using heudiconv.

## Features

- Query DICOM studies from CFMM with flexible search specifications
- Filter studies with include/exclude rules
- Download DICOM studies from CFMM
- Convert DICOM to BIDS format using heudiconv
- Apply post-conversion fixes (remove files, update JSON metadata, fix NIfTI orientation)
- Validate BIDS datasets
- Generate quality control (QC) reports for each subject/session

## Workflow Stages

The workflow is organized into 5 main processing stages plus a final copy stage, each producing intermediate outputs:

### 1. Query Stage (`results/0_query`)
Queries DICOM studies from CFMM using search specifications defined in `config/config.yaml`. Features include:
- Multiple search specifications with different query parameters
- Flexible metadata mapping (e.g., extract subject/session from PatientID, StudyDate)
- Pattern matching with regex extraction
- Automatic sanitization of subject/session IDs
- Validation of subject/session ID format (alphanumeric only)

Output: `studies.tsv` - Complete list of matched studies

### 2. Filter Stage (`results/1_filter`)
Post-filters the queried studies based on include/exclude rules. Features include:
- Include/exclude filters using pandas query syntax
- Optional `--config head=N` to process only first N subjects for testing

Output: `studies_filtered.tsv` - Filtered list of studies to process

### 3. Download Stage (`results/2_download`)
Downloads DICOM studies from CFMM using `cfmm2tar`.

Output: `dicoms/sub-*/ses-*/` - Downloaded DICOM files

### 4. Convert Stage (`results/3_convert`)
Converts DICOMs to BIDS format using heudiconv and generates QC reports. Features include:
- BIDS conversion with heudiconv and custom heuristic
- QC report generation (series list and unmapped summary)
- BIDS validation with `bids-validator-deno`
- Metadata preservation (auto.txt, dicominfo.tsv, filegroup.json)

Outputs:
- `bids/sub-*/ses-*/` - BIDS-formatted data
- `info/sub-*/ses-*/` - Heudiconv metadata files
- `qc/sub-*/ses-*/` - QC reports (series.svg, unmapped.svg)
- `qc/bids_validator.json` - BIDS validation results

### 5. Fix Stage (`results/4_fix`)
Applies post-conversion fixes to the BIDS dataset. Available fix actions:
- **remove**: Remove files matching a pattern (e.g., unwanted fieldmaps)
- **update_json**: Update JSON sidecar metadata (e.g., add PhaseEncodingDirection)
- **fix_orientation**: Reorient NIfTI files to canonical RAS+ orientation

Outputs:
- `bids/sub-*/ses-*/` - Fixed BIDS data
- `info/sub-*/ses-*/sub-*_ses-*_provenance.json` - Fix provenance tracking
- `qc/bids_validator.json` - Post-fix BIDS validation results

### 6. Final Stage (`bids/`)
Copies the validated and fixed BIDS dataset to the final output directory.

## QC Reports

The workflow automatically generates QC reports for each subject/session after heudiconv conversion. The reports include:

1. **Series List** (`*_series.svg`): A detailed table showing each series with:
   - Series ID and description
   - Protocol name
   - Image dimensions
   - TR and TE values
   - Corresponding BIDS filename (or "NOT MAPPED" if unmapped)

2. **Unmapped Summary** (`*_unmapped.svg`): A summary of series that were not mapped to BIDS, helping identify potential missing data or heuristic issues

QC reports are saved in the convert stage: `results/3_convert/qc/sub-{subject}/ses-{session}/`

**Note:** The QC report generation is integrated into the Snakemake workflow as a script directive and cannot be run manually as a standalone CLI tool.

### Aggregate HTML Report

After the fix stage, an aggregate HTML report is automatically generated that consolidates QC information from all subjects and sessions. The report includes:

- **Overview Statistics**: Total subjects, sessions, series, and unmapped series count
- **BIDS Validation Results**: Validation results from both convert and fix stages
- **Aggregated Series Table**: All series data sorted by subject and session
- **Post-Conversion Fix Provenance**: Details of fixes applied to each session
- **Heudiconv Filegroup Metadata**: Detailed metadata from heudiconv conversion (collapsible)

The aggregate report is located at: `results/4_fix/qc/aggregate_report.html`

This report provides a comprehensive overview of the entire dataset conversion process and is useful for quality control and troubleshooting.

## Configuration

The workflow is configured via `config/config.yaml`. Key configuration sections include:

For working examples, see:
- `config/config_trident15T.yml` - Configuration for Trident 15T scanner
- `config/config_cogms.yml` - Configuration for CogMS study

### Query Configuration (`search_specs`)
Define one or more DICOM queries with metadata mappings:
```yaml
search_specs:
  - dicom_query:
      study_description: YourStudy^*
      study_date: 20230101-
    metadata_mappings:
      subject:
        source: PatientID        # Extract from PatientID field
        pattern: '_([^_]+)$'    # Regex to extract subject ID
        sanitize: true          # Remove non-alphanumeric characters
      session:
        source: StudyDate       # Use StudyDate as session ID
```

### Filter Configuration (`study_filter_specs`)
Post-filter studies with include/exclude rules:
```yaml
study_filter_specs:
  include:
    - "subject.str.startswith('sub')"  # Include only subjects starting with 'sub'
  exclude:
    - "StudyInstanceUID == '1.2.3.4.5'"  # Exclude specific study
```

### Download Configuration
- `cfmm2tar_download_options`: Options passed to cfmm2tar (e.g., `--skip-derived`)
- `credentials_file`: Path to CFMM credentials file
- `merge_duplicate_studies`: If `true`, automatically merge multiple studies for the same subject/session (default: `false`)

#### Merging Multiple Studies

When `merge_duplicate_studies: true` is enabled and multiple studies match the same subject/session:
- Each study's DICOM tar file is processed separately with heudiconv
- Outputs are automatically merged into a single session
- Series IDs are offset (incremented by 1000) to prevent conflicts between studies
- A `study_uid` column is added to track which series came from which study
- This is useful when subjects have multiple scan sessions on the same day (e.g., due to console reboot)

### Convert Configuration
- `heuristic`: Path to heudiconv heuristic file
- `dcmconfig_json`: Path to dcm2niix configuration
- `heudiconv_options`: Additional heudiconv options

#### Available Heuristics

The workflow includes several heuristic files for different scanner configurations:

- **`heuristics/cfmm_base.py`**: Base CFMM heuristic supporting standard sequences including:
  - MP2RAGE, MEMP2RAGE, Sa2RAGE
  - T2 TSE (Turbo Spin Echo)
  - T2 SPACE, T2 FLAIR
  - Multi-echo GRE
  - TOF Angiography
  - Diffusion-weighted imaging
  - BOLD fMRI (multiband, psf-dico)
  - Field mapping (EPI-PA, GRE)
  - **DIS2D/DIS3D distortion-corrected reconstructions** - Improved detection that robustly identifies distortion-corrected images regardless of their position in the DICOM image_type metadata
  
- **`heuristics/trident_15T.py`**: Trident 15T scanner-specific heuristic
- **`heuristics/Menon_CogMSv2.py`**: CogMS study-specific heuristic

The heuristics automatically detect and label distortion-corrected (DIS2D/DIS3D) reconstructions using the `rec-DIS2D` or `rec-DIS3D` BIDS suffix.

### Fix Configuration (`post_convert_fixes`)
Define fixes to apply after conversion:
```yaml
post_convert_fixes:
  - name: remove_fieldmaps
    pattern: "fmap/*dir-AP*"
    action: remove
    
  - name: add_phase_encoding
    pattern: "func/*bold.json"
    action: update_json
    updates:
      PhaseEncodingDirection: "j-"
      
  - name: reorient_nifti
    pattern: "anat/*T1w.nii.gz"
    action: fix_orientation
```

### Other Options
- `final_bids_dir`: Final output directory (default: `bids`)
- `stages`: Customize intermediate stage directories

## Usage

1. Install [pixi](https://pixi.sh/latest/installation/)
   ```bash
   curl -fsSL https://pixi.sh/install.sh | sh
   ```
2. Clone the [cfmm2bids repository](https://github.com/akhanf/cfmm2bids)
   ```bash
   git clone https://github.com/akhanf/cfmm2bids
   cd cfmm2bids
   ```
3. Install dependencies into pixi virtual environment
   ```bash
   pixi install
   ```
4. Configure your search specifications by editing the `config/config.yaml`
   
   **Note:** Example configurations are available:
   - `config/config_trident15T.yml` - Trident 15T scanner setup
   - `config/config_cogms.yml` - CogMS study setup
   
   You can use these as starting points or use the template in `config/config.yaml`.
   
   To use one of the example configs directly:
   ```bash
   pixi run snakemake --configfile config/config_trident15T.yml --dry-run
   ```

5. Run the workflow as a dry-run to see what will be executed:
   ```bash
   pixi run snakemake --dry-run
   ```

6. Run specific workflow stages or the full workflow:
   
   **Run all stages (query → filter → download → convert → fix → final)**:
   ```bash
   pixi run snakemake --cores all
   ```
   
   **Run only download stage**:
   ```bash
   pixi run snakemake download --cores all
   ```
   
   **Run only convert stage** (includes QC reports):
   ```bash
   pixi run snakemake convert --cores all
   ```
   
   **Run only fix stage**:
   ```bash
   pixi run snakemake fix --cores all
   ```
   
   **Process only the first subject** (useful for testing):
   ```bash
   pixi run snakemake head --cores all
   ```
   
   **Process only first N subjects** (e.g., first 3):
   ```bash
   pixi run snakemake --config head=3 --cores all
   ```

7. Run the workflow on a SLURM cluster:
   ```bash
   pixi run snakemake --executor slurm --jobs 10
   ```   

## Output Directory Structure

### Final BIDS Output
```
bids/                           # Final BIDS-formatted output
├── dataset_description.json    # BIDS dataset metadata
└── sub-*/
    └── ses-*/                  # Subject/session data (anat/, func/, fmap/, etc.)
```

### Intermediate Stages
```
results/
├── 0_query/
│   └── studies.tsv                    # All queried studies
├── 1_filter/
│   └── studies_filtered.tsv           # Filtered studies to process
├── 2_download/
│   └── dicoms/
│       └── sub-*/ses-*/               # Downloaded DICOM files
├── 3_convert/
│   ├── bids/
│   │   └── sub-*/ses-*/               # Initial BIDS conversion
│   ├── info/
│   │   └── sub-*/ses-*/               # Heudiconv metadata (auto.txt, dicominfo.tsv, etc.)
│   └── qc/
│       ├── sub-*/ses-*/               # QC reports (series.svg, unmapped.svg)
│       └── bids_validator.json        # BIDS validation results
└── 4_fix/
    ├── bids/
    │   └── sub-*/ses-*/               # Fixed BIDS data
    ├── info/
    │   └── sub-*/ses-*/               # Fix provenance files
    └── qc/
        ├── bids_validator.json        # Post-fix validation results
        └── aggregate_report.html      # Aggregate QC report for all sessions
```

## Repository Directory Structure
```
├── workflow/                   # Workflow files
│   ├── Snakefile              # Main Snakemake workflow
│   ├── lib/                   # Python modules
│   │   ├── query_filter.py   # DICOM query and filtering functions
│   │   ├── bids_fixes.py     # Post-conversion fix implementations
│   │   └── utils.py          # Utility functions
│   └── scripts/               # Workflow scripts
│       ├── generate_convert_qc_figs.py  # QC report generation
│       └── post_convert_fix.py          # Post-conversion fix application
├── heuristics/                 # Heudiconv heuristic files
│   ├── cfmm_base.py           # Base CFMM heuristic (supports DIS2D/DIS3D reconstruction)
│   ├── trident_15T.py         # Trident 15T scanner-specific heuristic
│   └── Menon_CogMSv2.py       # CogMS study-specific heuristic
├── resources/                  # Resource files
│   ├── dcm2niix_config.json   # dcm2niix configuration
│   └── dataset_description.json  # BIDS dataset metadata template
├── config/                     # Configuration files
│   ├── config.yaml            # Configuration template (customize this)
│   ├── config_trident15T.yml  # Example: Trident 15T scanner configuration
│   └── config_cogms.yml       # Example: CogMS study configuration
└── pixi.toml                  # Pixi project configuration and dependencies
```

## Available Target Rules

The workflow provides several target rules for running specific stages:

- `all` (default): Run all stages from query to final BIDS output
- `head`: Process only the first subject (useful for testing)
- `download`: Run only query, filter, and download stages
- `convert`: Run through convert stage (includes download, conversion, and QC)
- `fix`: Run through fix stage (includes all above plus post-conversion fixes)

Example usage:
```bash
# Test with first subject only
pixi run snakemake head --cores all

# Download all DICOMs without conversion
pixi run snakemake download --cores all

# Convert and generate QC reports
pixi run snakemake convert --cores all
```

