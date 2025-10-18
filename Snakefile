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
targets = expand(
    'bids/sub-{subject}/ses-{session}',
    zip,
    subject=df.subject,
    session=df.session
)


rule all:
    input: targets


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
    output:
        bids_subj_dir=directory('bids/sub-{subject}/ses-{session}')
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
        )



