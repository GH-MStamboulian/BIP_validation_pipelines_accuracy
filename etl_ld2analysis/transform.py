"""
Prepare the line data for data analysis.
"""
import numpy as np
import pandas as pd


EXTRA_CUTOFFS = {'xtr_quant': [5.0, np.Inf], 'en_quant': [0.35, 130]}


def transform(extracted_data, use_rrr):
    """
    Take the in-memory line data and transform it to a formaat that is easier to analyze.

    Args:
        extracted_data (dict): the values are data frames, the data
        use_rrr (bool): if ``True``, then only look at the reportable range

    Returns:
        dict: the transformed data
        * ``manifest`` (pandas.DataFrame):
        * ``sample_qc`` (pandas.DataFrame):
        * ``calls`` (pandas.DataFrame):
    """
    col_mapper = {
        'Patient ID': 'patient_id',
        'Sample ID': 'run_sample_id',
        'Flowcell ID': 'runid',
        'Total Molecules': 'fams',
        'Nucleotide change': 'mut_nt',
        'Alteration Type': 'variant_type',
        'Clinical Variant Class': 'clinical_class',
        'Unique Variant Identifier': 'variant_key_nt',
        'MAF/CN': 'percentage',
        'Mutant Molecules': 'mut_fams',
        'Call': 'call',
        'Amino acid change': 'mut_aa',
        'Genomic position': 'position',
        'Gene': 'gene',
        'Sequencing Data Folder': 'data_location',
        'Raw Sequencing Data Folder': 'raw_data_location',
        'Processing Lab': 'assay',
        'Panel': 'panel',
        'Biobank Source': 'biobank',
        'Sample Collection': 'collection',
        'Detection': 'pa_call',
        'Extraction Yield': 'xtr_quant',
        'Enrichment Molarity': 'en_quant',
        'QC Note': 'exclude_reason',
        'Cancer Type': 'cancer_type',
        'In RRR': 'in_rrr'
    }
    df_manifest = extracted_data['manifest']
    df_manifest = df_manifest.rename(columns=col_mapper)

    # The QC data
    df_qc = extracted_data['sample_qc']
    qc_table = df_qc.rename(columns=col_mapper)
    df_cdx = df_manifest[df_manifest['assay'] == 'cdx']
    sample_qc_table = qc_table[qc_table['run_sample_id'].isin(df_cdx['run_sample_id'])]

    index_columns = ['runid', 'run_sample_id']
    extra_columns = ['xtr_quant', 'en_quant']
    df_extra = df_manifest[['patient_id', 'assay', 'xtr_quant', 'en_quant']]
    for extra_column in extra_columns:
        low_cutoff, high_cutoff = EXTRA_CUTOFFS[extra_column]
        df_extra[f'{extra_column}_status'] = \
            ((df_extra[extra_column] > low_cutoff) & (df_extra[extra_column] < high_cutoff)).\
            replace({True: 'PASS', False: 'FAIL'})
    df_extra[extra_columns] = df_extra[extra_columns].round(2)
    # Make the tall QC table a wide table
    df_extra = pd.merge(df_extra, df_cdx[['runid', 'run_sample_id',
                                          'patient_id', 'assay']], how='left')
    tdf_sample_qc = df_qc_to_ld(sample_qc_table, 'sample', dict(), index_columns, df_extra)
    # Add the "exclude_reason" as a pass/fail metric
    df_manifest['exclude_reason'] = df_manifest['exclude_reason'].fillna('').astype(str)
    # Right merge with the manifest and fill in total QC failure, pass/fail for non-CDx samples
    tdf_sample_qc = pd.merge(
        tdf_sample_qc,
        df_manifest[['runid', 'run_sample_id', 'patient_id', 'assay', 'exclude_reason']],
        how='right')
    tdf_sample_qc['exclude_reason_status'] = \
        tdf_sample_qc['exclude_reason'].apply(lambda x: 'PASS' if x == '' else 'FAIL')
    status_columns = [col for col in tdf_sample_qc if col.endswith('_status')]
    tdf_sample_qc[status_columns] = tdf_sample_qc[status_columns].fillna('NODATA')

    # Make a single "total" column that records if all metrics passed QC or if any failed
    # The germline contamination metric is actually a review state
    # https://guardanthealth.atlassian.net/wiki/spaces/GIVDTECH/pages/878478780
    special_status = ['sample_autoqc_total_status', 'sample_germline_contamination_status',
                      'exclude_reason_status']
    extra_status = [col + '_status' for col in extra_columns]
    regular_status = [col for col in status_columns
                      if col not in extra_status + special_status]
    passed_regular = (tdf_sample_qc[regular_status] != 'FAIL').all(axis=1)
    passed_extra = (tdf_sample_qc[extra_status] != 'FAIL').all(axis=1)
    passed_bip = tdf_sample_qc['sample_autoqc_total_status'] != 'FAIL'
    passed_contam = tdf_sample_qc['sample_germline_contamination_status'] != 'FAIL'
    passed_exclude = tdf_sample_qc['exclude_reason_status'] != 'FAIL'
    passed_qc = passed_exclude & passed_extra & (passed_bip | (passed_regular & ~passed_contam))
    tdf_sample_qc['total'] = passed_qc.replace({True: 'PASS', False: 'FAIL'})
    tdf_sample_qc = pd.merge(tdf_sample_qc,
                             df_manifest[['patient_id', 'assay', 'collection']], how='left')

    # The call data
    df_calls = extracted_data['calls']
    df_calls = df_calls.rename(columns=col_mapper)
    df_calls['clinical_class'] = df_calls['clinical_class'].fillna('')
    df_calls['call'] = df_calls['call'].replace({'None': None})

    # Label if something is a clinically relevant category
    df_calls['is_cs'] = df_calls['clinical_class'] != ''

    # Annotate the calls with the manifeest information
    df_calls = pd.merge(df_calls, df_manifest, how='left')

    # Only look at varianats in the RRR
    if use_rrr:
        df_calls = df_calls[df_calls['in_rrr'].values]

    # The reviewed call data
    df_reviewed = extracted_data['reviewed_calls'].rename(columns=col_mapper)

    # The CNV call data
    df_cnv_calls = extracted_data['cnv_calls']
    df_cnv_calls = df_cnv_calls.rename(columns=col_mapper)
    df_cnv_calls = pd.merge(df_cnv_calls, df_manifest, how='left')

    out = {
        'manifest': df_manifest,
        'sample_qc': tdf_sample_qc,
        'calls': df_calls,
        'reviewed_calls': df_reviewed,
        'cnv_calls': df_cnv_calls
    }
    return out


def df_qc_to_ld(qc_table, category, column_mapper, index_columns, df_extra=None):
    """
    Convert the sample QC table (a "tall" data frame) to a "wide" data frame where each column has
    the data values for a QC metric. The names of the metric are given a "_status" suffix to record
    the status -- PASS/FAIL/REVIEW.
    Args:
        qc_table (pandas.DataFrame): QC tall table (like in GHDB)
        category (str): flowcell, sample, or control
        column_mapper (dict): how to rename columns
        index_columns (list): the columns to use as the index
        df_extra (pandas.DataFrame): any extra QC information to merge in. It should already be in
            wide format. This could be a manifest of samples so that samples without any QC data
            will have null values when making the final table.
    Returns:
        pandas.DataFrame: the data ready for loading into a tab
    """
    # Reshape from tall to wide data frame
    df_category = qc_table[qc_table['category'] == category]
    df_values = pd.pivot_table(
        df_category, index=index_columns, columns='metric', values='value', aggfunc=lambda x: x)
    df_status = pd.pivot_table(
        df_category, index=index_columns, columns='metric', values='status', aggfunc=lambda x: x)
    assert df_status.shape == df_values.shape

    # Combine the PASS/FAIL status with the values by appending suffix to status columns
    df_status.columns = [col + '_status' for col in df_status.columns]
    df_qc = pd.concat([df_values, df_status], axis=1).reset_index()
    assert len(df_qc) == len(df_status)

    # Add in extra columns if provided
    if df_extra is not None:
        df_qc = pd.merge(df_qc, df_extra, how='outer', on=index_columns)
        df_qc = df_qc.fillna('NA')

    # Rename columns and then done
    df_qc = df_qc.rename(columns=column_mapper)
    return df_qc
