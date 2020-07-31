"""
The protocol calls for us to to calculate NPA for different classes of variants.
"""
import os
import pandas as pd
from joblib import Memory
import yaml
from av360.extract.extract_calls import get_catted_call_tables
from av360.historical.identifiers import get_past_unique_ids
from settings import CACHE_DIR, DATA_DIR
from etl_ld2analysis.load import load as load_study_data
from etl_gh2ld.transform_calls import annotate_homopolymers, annotate_long_indels


memory = Memory(CACHE_DIR, verbose=0)
c_get_catted_call_tables = memory.cache(get_catted_call_tables)
c_get_past_unique_ids = memory.cache(get_past_unique_ids)

# ------------------ Read in the cleaned up data --------------------------------------------------
use_rrr = True
if use_rrr:
    study_data = load_study_data(
        line_data_path=os.path.join(DATA_DIR, 'line_data_new.xlsx'), use_rrr=True)
else:
    study_data = load_study_data(
        line_data_path=os.path.join(DATA_DIR, 'line_data_old.xlsx'), use_rrr=False)
df_calls = study_data['calls']

# ------------------- Get previous calls ----------------------------------------------------------
df_past = c_get_past_unique_ids(
    'ghdb_prod', start_date='2018-01-01', end_date='2018-12-31', panels=['GH2.10.1'],
    bip_versions=['3.3'], production=True, qc_pass=True, patient_info=True)
pre_filters = {'pos_only': True,
               'pos_cnv_labels': tuple([1, 2, 3]),
               'has_rsid_suffix': True}
post_filters = {
    'av_variants_only': True,
    'reportable_only': True,
    'collapse_fusions': True
}
df_past_calls = c_get_catted_call_tables(df_past['unique_id'], 'ghdb_prod',
                                         pre_filters=pre_filters,
                                         post_filters=post_filters)
df_past_calls['homopolymer'] = annotate_homopolymers(df_past_calls)
df_past_calls['long_indel'] = annotate_long_indels(df_past_calls)


# ------------ Number of variants evaluates for NPA for each variant category ---------------------
n_panel = df_calls[['variant_type', 'variant_key_nt']].drop_duplicates().\
    groupby('variant_type').size()
df_cs = df_calls[df_calls['is_cs']]
n_cs = df_cs[['variant_type', 'variant_key_nt']].drop_duplicates().\
    groupby('variant_type').size()

N_VAR_SPACE = pd.Series({(True, 'snv'): n_cs['snv'],
                         (False, 'snv'): n_panel['snv'],
                         (True, 'indel'): n_cs['indel'],
                         (False, 'indel'): n_panel['indel'],
                         (True, 'cnv'): 2,
                         (True, 'fusion'): 4})
N_VAR_SPACE.index.names = ['is_cs', 'variant_type']
N_VAR_SPACE.name = 'n_var'

# Also save what is hard-coded into Table 10 of D-000111
N_VAR_SPACE_PROTOCOL = pd.Series({
    (True, 'snv'): 1166,
    (False, 'snv'): 36744,
    (True, 'indel'): 685,
    (False, 'indel'): 8101,
    (True, 'cnv'): 2,
    (True, 'fusion'): 4})
N_VAR_SPACE_PROTOCOL.index.names = ['is_cs', 'variant_type']
N_VAR_SPACE_PROTOCOL.name = 'n_var_protocol'

# ----------------- Calculate priors --------------------------------------------------------------
priors = dict()

# First, clinically significant variants
n_patients = len(df_past)
df_counts = df_past_calls.groupby('clinical_class').\
    apply(lambda x: pd.Series({'n_var': len(set(x['variant_key_nt'])),
                               'n_pos': len(x)}))
df_counts.loc[['alk_fusion', 'ntrk1_fusion', 'ret_fusion', 'ros1_fusion',
               'egfr_l858r', 'egfr_t790m'], 'n_var'] = 1
df_counts['prior'] = df_counts['n_pos'] / (df_counts['n_var'] * n_patients)
df_counts['var_per_sample'] = df_counts['n_pos'] / n_patients
df_counts = df_counts.drop('')
priors['clinical_class'] = dict()
priors['clinical_class']['n_var'] = df_counts['n_var'].to_dict()
priors['clinical_class']['prior'] = df_counts['prior'].to_dict()
priors['clinical_class']['var_per_sample'] = df_counts['var_per_sample'].to_dict()


# Then, long indels
df_counts = df_past_calls.groupby('long_indel').\
    apply(lambda x: pd.Series({'n_var': len(set(x['variant_key_nt'])),
                               'n_pos': len(x)}))
df_counts['prior'] = df_counts['n_pos'] / (df_counts['n_var'] * n_patients)
df_counts['var_per_sample'] = df_counts['n_pos'] / n_patients
df_long_indel = df_counts.drop(False)
df_long_indel.index = ['long_indel']

# Then, homopolymers
df_counts = df_past_calls.groupby('homopolymer').\
    apply(lambda x: pd.Series({'n_var': len(set(x['variant_key_nt'])),
                               'n_pos': len(x)}))
df_counts['prior'] = df_counts['n_pos'] / (df_counts['n_var'] * n_patients)
df_counts['var_per_sample'] = df_counts['n_pos'] / n_patients
df_homopolymer = df_counts.drop(False)
df_homopolymer.index = ['homopolymer']
df_other = pd.concat([df_long_indel, df_homopolymer])
priors['special'] = dict()
priors['special']['n_var'] = df_other['n_var'].to_dict()
priors['special']['prior'] = df_other['prior'].to_dict()
priors['special']['var_per_sample'] = df_other['var_per_sample'].to_dict()

# ------------------ Append priors from protocol --------------------------------------------------
#: Number of calls per variant from historical data
LDT_PRIORS = pd.Series({(True, 'snv'): 0.000232,
                        (False, 'snv'): 0.000059,
                        (True, 'indel'): 0.000089,
                        (False, 'indel'): 0.000040,
                        (True, 'cnv'): 0.116904,
                        (True, 'fusion'): 0.002542})
LDT_PRIORS.index.names = ['is_cs', 'variant_type']
LDT_PRIORS.name = 'prior'

FREQ_PRIORS = pd.Series({(True, 'snv'): 0.271,
                         (False, 'snv'): 2.182,
                         (True, 'indel'): 0.061,
                         (False, 'indel'): 0.322,
                         (True, 'cnv'): 0.234,
                         (True, 'fusion'): 0.010})
FREQ_PRIORS.index.names = ['is_cs', 'variant_type']
FREQ_PRIORS.name = 'var_per_sample'

df_primary = pd.concat([N_VAR_SPACE, N_VAR_SPACE_PROTOCOL, LDT_PRIORS, FREQ_PRIORS], axis=1).reset_index()
df_primary['type'] = 'variant_category'
df_primary['group'] = list(zip(df_primary['is_cs'], df_primary['variant_type']))
df_primary = df_primary.set_index('group')
priors['category'] = dict()
priors['category']['n_var'] = df_primary['n_var'].to_dict()
priors['category']['n_var_protocol'] = df_primary['n_var_protocol'].to_dict()
priors['category']['prior'] = df_primary['prior'].to_dict()
priors['category']['var_per_sample'] = df_primary['var_per_sample'].to_dict()

# -------- Number of variants evaluated for NPA for each variant class ----------------------------
# Update to use just the accuracy calls and not the numbers from historical data or from the
# protocol.
n_long_indel = df_calls[df_calls['long_indel']]['variant_key_nt'].nunique()
n_homopolymer = df_calls[df_calls['homopolymer']]['variant_key_nt'].nunique()
priors['special']['n_var'] = {'long_indel': n_long_indel, 'homopolymer': n_homopolymer}

df_cs = df_calls[['variant_key_nt', 'clinical_class', 'variant_key_nt']].drop_duplicates()
n_cs = df_cs['clinical_class'].value_counts()
single_var_classes = ['alk_fusion', 'ntrk1_fusion', 'ret1_fusion', 'ros_fusion',
                      'egfr_l858r', 'egfr_t790m', 'met_amplification', 'erbb2_amplification']
n_cs_update = {key: value for (key, value) in n_cs.items()
               if key not in single_var_classes and key != ''}
priors['clinical_class']['n_var'].update(n_cs_update)

# ---------- Save all priors in a single settings file --------------------------------------------
suffix = '_rrr' if use_rrr else ''
with open(os.path.join(DATA_DIR, f'priors{suffix}.yaml'), 'w') as file_pointer:
    print(yaml.dump(priors), file=file_pointer)
