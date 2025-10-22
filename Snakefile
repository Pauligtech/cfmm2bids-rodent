import pandas as pd
from pathlib import Path
from cfmm2tar import query_metadata, download_studies

configfile: 'config.yml'

all_dfs = []

for spec in config['search_specs']:
    # Query DICOM metadata
    df_ = query_metadata(
        credentials_file=config['credentials_file'],
        return_type="dataframe",
        **spec['dicom_query']
    )

    # Skip if empty
    if df_.empty:
        continue

    # Apply metadata extraction settings
    mappings = spec.get('metadata_mappings', {})
    for target, mapping in mappings.items():
        source_col = mapping['source']
        series = df_[source_col]

        # Optional regex extraction
        if 'pattern' in mapping:
            series = series.str.extract(mapping['pattern'], expand=False)

        # Optional cleaning / sanitization
        if mapping.get('sanitize', True):
            series = series.str.replace(r'[^A-Za-z0-9]', '', regex=True)

        # Optional remapping of specific values
        if 'map' in mapping:
            series = series.replace(mapping['map'])

        # Assign to target field
        df_[target] = series

    # Record query info for traceability (optional)
    df_['query_params'] = str(spec['dicom_query'])

    all_dfs.append(df_)

# Combine all query results into a single DataFrame
df = pd.concat(all_dfs, ignore_index=True)

localrules: dataset_description

# Build BIDS-style output targets
sub_ses_targets = expand(
    'bids/sub-{subject}/ses-{session}',
    zip,
    subject=df.subject,
    session=df.session
)

# Build QC report targets
qc_report_targets = expand(
    'qc/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_gantt.svg',
    zip,
    subject=df.subject,
    session=df.session
)


rule all:
    input: 
        sub_ses_targets,
        'bids/dataset_description.json',
        qc_report_targets


def get_uid_from_wildcards(wildcards):
    return lookup(
        query=f"subject == '{wildcards.subject}' and session == '{wildcards.session}'",
        within=df
    ).StudyInstanceUID



rule download_tar:
    params:
        uid=get_uid_from_wildcards
    output:
        dicoms_dir=directory('sourcedata/sub-{subject}/ses-{session}')
    log:
        'logs/download_tar/sub-{subject}_ses-{session}.log'
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15
    run:
        import sys
        from pathlib import Path
        
        # Ensure log directory exists
        Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        
        # Redirect stdout and stderr to log file
        with open(log[0], 'w') as log_file:
            original_stdout = sys.stdout
            original_stderr = sys.stderr
            try:
                sys.stdout = log_file
                sys.stderr = log_file
                download_studies(
                    output_dir=output.dicoms_dir,
                    credentials_file=config['credentials_file'],
                    study_instance_uid=params.uid,
                    **config['download_options']
                )
            finally:
                sys.stdout = original_stdout
                sys.stderr = original_stderr


rule heudiconv:
    input:
        dicoms_dir='sourcedata/sub-{subject}/ses-{session}',
        heuristic=config['heuristic'],
        dcmconfig_json=config['dcmconfig_json'],
    params: 
        heudiconv_options=config['heudiconv_options'],
        in_auto_txt='bids/.heudiconv/{subject}/ses-{session}/info/{subject}_ses-{session}.auto.txt',
        in_dicominfo_tsv='bids/.heudiconv/{subject}/ses-{session}/info/dicominfo_ses-{session}.tsv',
        in_filegroup_json='bids/.heudiconv/{subject}/ses-{session}/info/filegroup_ses-{session}.json',
        out_info_dir='sourcedata/heudiconv/sub-{subject}/ses-{session}',
    output:
        bids_subj_dir=directory('bids/sub-{subject}/ses-{session}'),
        auto_txt='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_auto.txt',
        dicominfo_tsv='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_dicominfo.tsv',
        filegroup_json='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_filegroup.json'
    log:
        'logs/heudiconv/sub-{subject}_ses-{session}.log'
    shadow: 'minimal'
    threads: 16
    resources: 
        mem_mb=8000,
        runtime=15
    group: 'convert'
    shell:
        (
            "mkdir -p $(dirname {log}) && "
            "heudiconv --files {input.dicoms_dir}"
            " -c dcm2niix"
            " -o bids"
            " -ss {wildcards.session}"
            " -s {wildcards.subject}"
            " -f {input.heuristic}"
            " --bids notop"
            " --dcmconfig {input.dcmconfig_json}"
            " --overwrite"
            " {params.heudiconv_options}"
            " && mkdir -p {params.out_info_dir}"
            " && cp {params.in_auto_txt} {output.auto_txt}"
            " && cp {params.in_dicominfo_tsv} {output.dicominfo_tsv}"
            " && cp {params.in_filegroup_json} {output.filegroup_json}"
            " &> {log}"
        )

rule dataset_description:
    input:
        'resources/dataset_description.json'
    output:
        'bids/dataset_description.json'
    log:
        'logs/dataset_description/dataset_description.log'
    shell:
        'mkdir -p $(dirname {log}) && cp {input} {output} &> {log}'


rule generate_qc_report:
    input:
        auto_txt='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_auto.txt',
        dicominfo_tsv='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_dicominfo.tsv',
        filegroup_json='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_filegroup.json'
    output:
        gantt='qc/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_gantt.svg',
        series_list='qc/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_series-list.svg',
        unmapped='qc/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_unmapped.svg'
    log:
        'logs/generate_qc_report/sub-{subject}_ses-{session}.log'
    threads: 1
    resources: 
        mem_mb=4000,
        runtime=10
    group: 'convert'
    script:
        'scripts/generate_qc_report.py'

