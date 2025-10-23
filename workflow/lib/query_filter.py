import pandas as pd
from cfmm2tar import query_metadata, download_studies

# Function to validate a column
def validate_column(df, col):
    # Non-blank: not null and not whitespace
    non_blank = df[col].notna() & df[col].str.strip().ne('')
    
    # Alphanumeric check
    pattern = r'^[A-Za-z0-9]+$'
    alphanumeric = df[col].str.match(pattern, na=False)
    
    # Combine both conditions
    valid = non_blank & alphanumeric
    
   
    return valid



def query_dicoms(search_specs, **query_metadata_kwargs):

    all_dfs = []

    for spec in search_specs:
        # Query DICOM metadata
        df_ = query_metadata(
            return_type="dataframe",
            **query_metadata_kwargs,
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

    return df

def post_filter(df, post_filter_specs):
    if not post_filter_specs:
        return df

    for q in (post_filter_specs.get("include") or []):
        df = df.query(q)

    for q in (post_filter_specs.get("exclude") or []):
        df = df.query(f"not ({q})")

    return df


