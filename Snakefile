import pandas as pd
from cfmm2tar import query_metadata, download_studies_from_metadata
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



#some processing (manual + regex + lookup etc) will produce the subject and session columns
#have done this manually in  study_metadata_withsubjses.tsv
df = pd.read_csv('study_metadata_withsubjses.tsv',sep='\t')
print(df)
print(df.subject)
targets = expand(bids(
             subject='{subject}',
             session='{session}',
             suffix='dicoms'), zip,
           subject=df.subject,
           session=df.session)

rule all:
    input: targets

rule download_tar:
    output:
        dicoms_dir=directory(bids(
             subject='{subject}',
             session='{session}',
             suffix='dicoms'))
    run:
        download_studies_from_metadata(
            output_dir=output.dicoms_dir,
            credentials_file='~/.uwo_credentials.bd',
            metadata=df_download,
            temp_dir='./temp_dicoms',
        )
