import pandas as pd
from cfmm2tar import query_metadata, download_studies
from snakebids import bids

if 0:
    # Get metadata as DataFrame for filtering
    df = query_metadata(
        credentials_file='~/.uwo_credentials.bd',
        study_description="Prado^*",
        study_date="20250910",
        return_type="dataframe"  # or "list" for list of dicts
    )
    #print(df)

    df_download = df.sort_values(by='StudyDate')
    df_download.to_csv('study_metadata.tsv',sep='\t',index=False)

#        download_studies_from_metadata(
#            output_dir=output.dicoms_dir,
#            credentials_file='~/.uwo_credentials.bd',
#            metadata=df_download,
#            temp_dir='./temp_dicoms',
#        )

#some processing (manual + regex + lookup etc) will produce the subject and session columns
#have done this manually in  study_metadata_withsubjses.tsv
df = pd.read_csv('study_metadata_withsubjses.tsv',sep='\t')
print(df)
print(df.subject)


targets = expand(bids(
             subject='{subject}',
             session='{session}',
             suffix='heudiconv.done'), zip,
           subject=df.subject,
           session=df.session)

rule all:
    input: targets

rule download_tar:
    params:
        uid=lambda wildcards: lookup(query=f"subject == '{wildcards.subject}'",
                   within=df).StudyInstanceUID
    output:
        dicoms_dir=directory(bids(
             subject='{subject}',
             session='{session}',
             suffix='dicoms'))
    run:
        print(params.uid)
        download_studies(
            output_dir=output.dicoms_dir,
            credentials_file='~/.uwo_credentials.bd',
            study_instance_uid=params.uid,
            temp_dir='./temp_dicoms',
        )

rule heudiconv:
    input:
        dicoms_dir=bids(
             subject='{subject}',
             session='{session}',
             suffix='dicoms'),
        heuristic='resources/heuristic.py',
        dcmconfig_json='resources/dcm2niix_config.json'

    output:
        flag=touch(bids(
             subject='{subject}',
             session='{session}',
             suffix='heudiconv.done'))

    shell:
        "heudiconv --files {input.dicoms_dir}"
        " -c dcm2niix"
        " -o ./test_bids"
        "  -ss {wildcards.session}"
        "  -s {wildcards.subject}"
        "  -f {input.heuristic}"
        "  --bids"
        "  --dcmconfig {input.dcmconfig_json} --overwrite"
        

