"""
Extract all of the data needed for the accuracy study and put it into memory.
"""
import os
import pandas as pd
from joblib import Memory
from D000111d.settings import DATA_DIR, CACHE_DIR, REPORT_DIR
from av360.utils.gsheets_connect import make_client
from av360.extract.extract_calls import get_catted_call_tables
from av360.extract.extract_qc_data import get_sample_qc_table
from av360.extract.extract_warehouse import get_acs_dump, extract_patient_info_table
from av360.extract.extract_flowcentral_files import read_globs_dataframe_remote
from av360.utils.wtf import get_raw_paths, get_processed_paths_info, get_processed_path
from av360.utils.utils import make_unique_ids
REVIEWED_CALLS_SHEET = '1wdl5JvDtDXOmDzsw5NidZVOkvJEGUDMmZHOgbm_7IaQ'


memory = Memory(CACHE_DIR, verbose=0)
c_get_catted_call_tables = memory.cache(get_catted_call_tables)
c_get_sample_qc_table = memory.cache(get_sample_qc_table)
c_get_acs_dump = memory.cache(get_acs_dump)
c_get_raw_paths = memory.cache(get_raw_paths)
c_get_processed_paths_info = memory.cache(get_processed_paths_info)
c_extract_patient_info_table = memory.cache(extract_patient_info_table)


def extract(use_old_bip=False):
    """
    Extract data from the various sources at Guardant and load them all into an object in memory.

    Returns:
        dict: just a container for data

        * ``manifest`` (pandas.DataFrame): the sample manifest
        * ``calls`` (pandas.DataFrame): concatenated call tables with additional columns
        * ``reviewed_calls`` (pandas.DataFrmae): information for manually adjudicated calls
        * ``qc`` (pandas.DataFrame): tall format QC data
        * ``quants`` (pandas.DataFrame): dump of ACS quant data for samples in this study
    """
    # Sample manifest
    df_manifest = pd.read_csv(os.path.join(DATA_DIR, 'sample_manifest.csv'))

    # Make sure that we are pulling calls from the correct version of BIP
    cdx_runids = df_manifest.loc[df_manifest['assay'] == 'cdx', 'runid'].unique()
    for runid in cdx_runids:
        flowcell_id = runid[-10:]
        df_paths = c_get_processed_paths_info(flowcell_id)
        if use_old_bip:
            df_old = df_paths[df_paths['bip_version'].str.contains('3.5.[01]')]
            df_old = df_old[df_old['size'] > 1000]
            assert not df_old.empty
            path = df_old.sort_values('date').tail(1)['path'].iloc[0]
            df_manifest.loc[df_manifest['runid'] == runid, 'data_location'] = path
        else:
            # runid's for samples run at GH
            df_cdx_runids = df_manifest.\
                loc[df_manifest['assay'] == 'cdx', ['runid', 'data_location']].\
                drop_duplicates()
            # Update the path to the BIP 3.5.3 data
            old_to_new_location = dict()
            for _, row in df_cdx_runids.iterrows():
                runid = row['runid']
                new_location = get_processed_path(
                    runid, bip_git_version='3.5.3-0-g8857b98', choose_biggest=True)
                assert isinstance(new_location, str)
                df_manifest.loc[df_manifest['runid'] == runid,
                                'bip_git_version'] = '3.5.3-0-g8857b98'
                df_manifest.loc[df_manifest['runid'] == runid,
                                'data_location'] = new_location
    df_manifest['unique_id'] = make_unique_ids(df_manifest['runid'], df_manifest['run_sample_id'],
                                               df_manifest['data_location'])

    # Calls
    # for the BIP 3.2 samples
    unique_ids_3p2 = \
        df_manifest.loc[df_manifest['bip_git_version'].str.startswith('3.2-0'), 'unique_id']
    # Reportability of regions being taken care of with bed files
    post_filters = {'reportable_only': False, 'av_variants_only': False, 'collapse_fusions': True}
    pre_filters = {'pos_only': False, 'has_rsid_suffix': True, 'is_bip3.2': True}
    print('getting BIP 3.2 calls')
    df_calls_3p2 = c_get_catted_call_tables(unique_ids_3p2, 'ssh', pre_filters=pre_filters,
                                            post_filters=post_filters)
    # and for everything BIP 3.3 or later versions
    print('getting BIP >3.2 calls')
    pre_filters.update({'is_bip3.2': False})
    unique_ids_n3p2 = df_manifest.loc[~df_manifest['unique_id'].isin(unique_ids_3p2), 'unique_id']
    df_calls_n3p2 = c_get_catted_call_tables(unique_ids_n3p2, 'ssh', pre_filters=pre_filters,
                                             post_filters=post_filters)
    # Combine everything together
    df_calls = pd.concat([df_calls_3p2, df_calls_n3p2], sort=False)

    # Get the reviewed calls
    gsheet_client = make_client()
    gsheet = gsheet_client.open_by_key(REVIEWED_CALLS_SHEET)
    sheet = gsheet.worksheet('final_adjustments')
    df_reviewed = pd.DataFrame(sheet.get_all_records())

    # QC data
    # Guardant QC data
    df_cdx = df_manifest[(df_manifest['assay'] == 'cdx') &
                         (df_manifest['exclude_reason'].isnull() |
                          (df_manifest['exclude_reason'] == ''))]
    df_runs = df_cdx[['data_location', 'bip_git_version']].drop_duplicates()
    sample_qc_table = c_get_sample_qc_table(df_runs['data_location'], 'sftp')

    # Get Guardant ACS data for as many samples as possible (the extraction yield data and the
    # enrichment quant)
    df_av = c_get_acs_dump('av')
    df_cv = c_get_acs_dump('cv')
    df_acs = pd.concat([df_av, df_cv]).drop_duplicates()
    df_quants = pd.merge(df_cdx[['run_sample_id']], 
                         df_acs[['run_sample_id', 'xtr_quant', 'en_quant']], how='left')
    df_quants = df_quants.dropna().drop_duplicates()

    # Get raw sequencing folders (study report must have them listed)
    df_cdx_ldt = df_manifest[df_manifest['assay'].isin(['cdx', 'ldt'])]
    raw_items = dict()
    for runid in df_cdx_ldt['runid'].dropna().unique():
        raw_paths = c_get_raw_paths(runid)
        if len(raw_paths) > 1:
            raw_paths = raw_paths.sort_values('size').tail(1)
        if len(raw_paths) == 0:
            # A hack since the WTF tool stopped working
            cdx_flowcells = df_manifest.loc[df_manifest['assay'] == 'cdx', 'runid']
            # assert runid in cdx_flowcells.values
            raw_paths = pd.DataFrame({'path': ['/ghds/ivd/raw/' + runid]})
        assert len(raw_paths) == 1
        raw_items[runid] = raw_paths['path'].iloc[0]
    df_raw_locations = pd.DataFrame(pd.Series(raw_items)).reset_index()
    df_raw_locations.columns = ['runid', 'raw_data_location']

    # Need to get cancer types
    old_dir = os.path.join(REPORT_DIR, '../D000111/data')
    coll1_path = os.path.join(old_dir, '100 Consecutive Samples List_08282018.xlsx')
    coll2_path = os.path.join(old_dir, 'Requested Variants Sample Details.xlsx')
    df1 = pd.read_excel(coll1_path).\
        rename(columns={'Cancer Type': 'cancer_type',
                        'BLIND #': 'run_sample_id'})
    df2 = pd.read_excel(coll2_path).\
        rename(columns={'Indication': 'cancer_type',
                        'BLIND #': 'run_sample_id'})
    df_patient = c_extract_patient_info_table()
    df_ldt = df_manifest[(df_manifest['collection'] == 3) & (df_manifest['assay'] == 'ldt')]
    df3 = pd.merge(df_ldt[['run_sample_id']], df_patient[['run_sample_id', 'cancer_type']],
                   how='left').rename(columns={'cancertype': 'cancer_type'})
    ct_columns = ['run_sample_id', 'cancer_type']
    df_ct = pd.concat([df1[ct_columns], df2[ct_columns], df3[ct_columns]])
    df_cancer_types = df_ct

    out = {
        'manifest': df_manifest,
        'calls': df_calls,
        'reviewed_calls': df_reviewed,
        'qc': sample_qc_table,
        'quants': df_quants,
        'raw_data_locations': df_raw_locations,
        'cancer_types': df_cancer_types
    }
    return out
