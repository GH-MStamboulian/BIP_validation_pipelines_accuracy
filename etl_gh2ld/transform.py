"""
Transform the extracted data into the format that will be saved as line data
"""
import pandas as pd
from av360.utils.utils import make_unique_ids
from D000111d.etl_gh2ld.transform_calls import transform_calls, transform_cnv_calls
from D000111d.etl_ld2analysis.transform import df_qc_to_ld, EXTRA_CUTOFFS


def transform(extracted_data):
    """
    Transform the extracted data from :py:func:`etl_gh2ld.extract` in preparation for excel
    spreadsheet.
    """
    # Sample Manifest
    column_mapper = {
        'patient_id': 'Patient ID',
        'run_sample_id': 'Sample ID',
        'runid': 'Flowcell ID',
        'data_location': 'Sequencing Data Folder',
        'panel': 'Panel',
        'assay': 'Processing Lab',
        'biobank': 'Biobank Source',
        'collection': 'Sample Collection',
        'call': 'Call',
        'pa_call': 'Detection',
        'variant_type': 'Alteration Type',
        'clinical_class': 'Clinical Variant Class',
        'gene': 'Gene',
        'position': 'Genomic position',
        'variant_key_nt': 'Unique Variant Identifier',
        'percentage': 'MAF/CN',
        'mut_fams': 'Mutant Molecules',
        'fams': 'Total Molecules',
        'mut_nt': 'Nucleotide change',
        'mut_aa': 'Amino acid change',
        'raw_data_location': 'Raw Sequencing Data Folder',
        'xtr_quant': 'Extraction Yield',
        'en_quant': 'Enrichment Molarity',
        'exclude_reason': 'QC Note',
        'cancer_type': 'Cancer Type',
        'in_rrr': 'In RRR',
        'rm_reportable': 'rm_reportable'
    }
    manifest_columns = ['patient_id', 'runid', 'run_sample_id', 'data_location',
                        'raw_data_location', 'panel', 'assay', 'biobank', 'collection',
                        'xtr_quant', 'en_quant', 'exclude_reason']
    df_manifest = extracted_data['manifest']

    # Remove samples that had an instrument failure at MDL (feedback from Trang)
    instrument_failure_patients = \
        df_manifest.loc[df_manifest['exclude_reason'] == 'Library Prep QC Fail', 'patient_id']
    df_manifest = df_manifest[~df_manifest['patient_id'].isin(instrument_failure_patients)]

    # Add the raw data location
    df_raw_locations = extracted_data['raw_data_locations']
    df_manifest = pd.merge(df_manifest, df_raw_locations, how='left')
    df_manifest['raw_data_location'] = df_manifest['raw_data_location'].fillna('')

    # Need to add the quant data here since it doesn't fit the format of QC line data
    df_quants = extracted_data['quants']
    df_manifest = pd.merge(df_manifest, df_quants, how='left')
    df_manifest[['xtr_quant', 'en_quant']] = \
        df_manifest[['xtr_quant', 'en_quant']].round(2).fillna('NA')
    df_manifest['exclude_reason'] = df_manifest['exclude_reason'].fillna('')
    tdf_manifest = df_manifest[manifest_columns].rename(columns=column_mapper)
    df_manifest['unique_id'] = make_unique_ids(df_manifest['runid'], df_manifest['run_sample_id'],
                                               df_manifest['data_location'])

    # Flowcell QC data
    qc_table = extracted_data['qc']
    tdf_flowcell_qc = qc_table[qc_table['category'] == 'flowcell'].copy()
    tdf_flowcell_qc = tdf_flowcell_qc.drop(['data_location', 'unique_id'], errors='ignore', axis=1)

    # The variant control QC
    tdf_control_qc = qc_table[qc_table['category'] == 'control'].copy()
    tdf_control_qc = tdf_control_qc.drop(['data_location', 'unique_id'], errors='ignore', axis=1)

    # The sample QC
    df_cdx = df_manifest[df_manifest['assay'] == 'cdx']
    tdf_sample_qc = qc_table[qc_table['category'] == 'sample'].copy()
    tdf_sample_qc = tdf_sample_qc[tdf_sample_qc['run_sample_id'].isin(df_cdx['run_sample_id'])]
    tdf_sample_qc = tdf_sample_qc.drop(['data_location', 'unique_id'], errors='ignore', axis=1)

    # Next, all of the calls
    df_calls = extracted_data['calls']
    df_reviewed = extracted_data['reviewed_calls']
    df_reviewed = df_reviewed[['patient_id', 'variant_type', 'variant_key_nt', 'cdx', 'mdl', 'ldt',
                               'cdx_reviewed', 'mdl_reviewed', 'ldt_reviewed']].copy()
    tdf_calls = transform_calls(df_calls, df_reviewed, df_manifest)

    # Create a data frame of the calls for other studies to use
    passed_run_sample_ids = find_passed_run_sample_ids(df_manifest, qc_table, extracted_data)
    df_dump = make_call_dump(tdf_calls, df_manifest, passed_run_sample_ids)

    # Next, all of the CNV calls -- *not* filtering out non-panel genes
    tdf_cnv_calls = transform_cnv_calls(df_calls)

    # Don't forget to rename columns for call information for line data
    tdf_calls = tdf_calls.rename(columns=column_mapper)
    tdf_cnv_calls = tdf_cnv_calls.rename(columns=column_mapper)

    # Fix up the cancer types and then merge it with the manifest
    df_ct = extracted_data['cancer_types']
    cancer_types = {
        'Bladder carcinoma': 'Bladder',
        'Breast': 'Breast',
        'Breast Cancer': 'Breast',
        'Carcinoma of unknown primary (CUP)': 'Other',
        'Cholangiocarcinoma': 'Colon',
        'Colorectal adenocarcinoma': 'Colon',
        'Cutaneous squamous cell carcinoma': 'Skin',
        'Endocrine': 'Endocrine',
        'Endocrine Tumor': 'Endocrine',
        'Esophageal/Gastroesophageal junction Adenocarcinoma': 'Esophagus',
        'GI': 'Gastrointestinal',
        'GI Tumor': 'Gastrointestinal',
        'GU': 'Genitourinary',
        'GU Tumor': 'Genitourinary',
        'GYN': 'Gynocelogic',
        'Gastric adenocarcinoma': 'Stomach',
        'Head and Neck': 'Head and Neck',
        'Head and neck/Thyroid': 'Head and Neck',
        'Hepatocellular carcinoma': 'Liver',
        'Large cell lung carcinoma': 'Lung',
        'Lung adenocarcinoma': 'Lung',
        'Lung cancer': 'Lung',
        'Lung squamous cell carcinoma': 'Lung',
        'Non-small cell lung carcinoma (NSCLC)': 'Lung',
        'Other': 'Other',
        'Ovarian carcinoma': 'Ovary',
        'Pancreatic ductal adenocarcinoma': 'Pancreas',
        'Prostate adenocarcinoma': 'Prostate',
        'Renal Cell Carcinoma': 'Kidney',
        'Renal pelvis urothelial carcinoma': 'Kidney',
        'Sarcoma': 'Other',
        'Small cell lung carcinoma': 'Lung',
        'Thoracic': 'Lung',
        'Thoracic Tumor': 'Lung'
    }
    df_ct['cancer_type'] = df_ct['cancer_type'].replace(cancer_types)
    df_ct = pd.merge(df_ct, df_manifest[['run_sample_id', 'patient_id']], how='left')
    df_ct = df_ct.drop('run_sample_id', axis=1)
    tdf_manifest = pd.merge(tdf_manifest, df_ct.rename(columns=column_mapper), how='left')

    transformed_data = {
        'flowcell_qc': tdf_flowcell_qc,
        'control_qc': tdf_control_qc,
        'sample_qc': tdf_sample_qc,
        'calls': tdf_calls,
        'manifest': tdf_manifest,
        'call_dump': df_dump,
        'reviewed_calls': df_reviewed.rename(columns=column_mapper),
        'cnv_calls': tdf_cnv_calls
    }
    return transformed_data


def make_call_dump(df_calls, df_manifest, passed_run_sample_ids):
    # Only include samples from specific assays
    included_assays = ['cdx', 'ldt', 'mdl']
    assay_run_sample_ids = df_manifest.loc[df_manifest['assay'].isin(included_assays),
                                           'run_sample_id']

    # Only include samples passing QC
    passed_runs = set(passed_run_sample_ids).intersection(assay_run_sample_ids)
    df_assay_calls = df_calls[df_calls['run_sample_id'].isin(passed_runs)]

    # Only include calls that are in specified panels
    panels = ['GH2.11', 'MDACC_v1.0', 'GH2.10.1']
    panel_columns = ['panel_' + panel for panel in panels]
    in_all_panels = df_assay_calls[panel_columns].all(axis=1)
    df_panel_calls = df_assay_calls[in_all_panels]

    # Only dump the positive calls
    df_dump = df_panel_calls[df_panel_calls['call'] > 0].copy()

    # Onluy the calls in the reportable range
    # df_dump = df_dump[df_dump['in_rrr']]
    return df_dump


def find_passed_run_sample_ids(df_manifest, qc_table, extracted_data):
    """
    I have to jump ahead a bit in the code to find which samples passed QC. This is some hacky
    copy and paste stuff.
    """
    # Sample QC data
    # Only look at samples in this study. There might be non-study samples on the flowcells
    df_cdx = df_manifest[df_manifest['assay'] == 'cdx']
    sample_qc_table = qc_table[qc_table['unique_id'].isin(df_cdx['unique_id'])]
    # Format the quant data, adding status columns
    index_columns = ['runid', 'run_sample_id']
    extra_columns = ['xtr_quant', 'en_quant']
    df_extra = extracted_data['quants']
    for extra_column in extra_columns:
        low_cutoff, high_cutoff = EXTRA_CUTOFFS[extra_column]
        df_extra[f'{extra_column}_status'] = \
            ((df_extra[extra_column] > low_cutoff) & (df_extra[extra_column] < high_cutoff)).\
            replace({True: 'PASS', False: 'FAIL'})
    df_extra[extra_columns] = df_extra[extra_columns].round(2)
    # Make the tall QC table a wide table
    df_extra = pd.merge(df_extra, df_cdx[['runid', 'run_sample_id']], how='left')
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
    # First line checks if theeir is any "force failure" because of a note of an instrument failure
    # from MDL or from us
    passed_qc = (passed_exclude &
                 # Second line checks if all of the quants passed for the CDx data
                 passed_extra &
                 # Third line checks of all of the BIP metrics passed *OR*
                 # if the germline contamination failed but all of the other BIP metrics passed
                 (passed_bip | (passed_regular & ~passed_contam)))
    tdf_sample_qc['total'] = passed_qc.replace({True: 'PASS', False: 'FAIL'})
    passed_run_sample_ids = tdf_sample_qc.loc[tdf_sample_qc['total'] == 'PASS', 'run_sample_id']
    return passed_run_sample_ids
