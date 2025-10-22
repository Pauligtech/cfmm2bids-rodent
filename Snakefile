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


localrules: download_tar

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
    run:
        download_studies(
            output_dir=output.dicoms_dir,
            credentials_file=config['credentials_file'],
            study_instance_uid=params.uid,
            **config['download_options']
        )


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
    shadow: 'minimal'
    threads: 16
    resources: 
        mem_mb=8000
    shell:
        (
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
        )

rule dataset_description:
    input:
        'resources/dataset_description.json'
    output:
        'bids/dataset_description.json'
    shell:
        'cp {input} {output}'


rule generate_qc_report:
    input:
        auto_txt='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_auto.txt',
        dicominfo_tsv='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_dicominfo.tsv',
        filegroup_json='sourcedata/heudiconv/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_filegroup.json'
    output:
        gantt='qc/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_gantt.svg',
        series_list='qc/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_series-list.svg',
        unmapped='qc/sub-{subject}/ses-{session}/sub-{subject}_ses-{session}_unmapped.svg'
    script:
        'scripts/generate_qc_report.py'

