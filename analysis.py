"""
The top level analysis script for the accuracy study.
"""
import sys
sys.path.append("/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/repositories/data_science/08_FDA/G360/180801_ANALYTICAL_VALIDATION/")
import os
from D000111d.etl_ld2analysis.load import load as load_study_data
from D000111d.compare_calls.compare_calls import get_call_comparisons, get_call_comparisons_2test
from D000111d.compare_calls.focal_cnvs import compare_focal_calls_coll123, compare_nonfocal_status
from D000111d.compare_calls.compare_groups import compare_collection_groups, compare_combined_collections
from D000111d.tables import table_qc_counts_v2, table_qc_details, table_call_count_summary, \
    tables_ppa_npa_collections_summary, tables_ppa_npa_combined_v4, \
    table_prior_and_theta_values, table_variant_override_info, table_cancer_types, \
    tables_ppa_npa_combined_bips
from D000111d.settings import DATA_DIR


line_file = sys.argv[1]
table_out_dir = sys.argv[2]

if (not os.path.isdir(table_out_dir)):
    os.makedirs(table_out_dir)

# ----- Load the study data -----------------------------------------------------------------------
#study_data = load_study_data(
#    line_data_path=os.path.join(DATA_DIR, 'line_data.xlsx'), use_rrr=False)
study_data = load_study_data(line_data_path=line_file, use_rrr = False)
df_manifest = study_data['manifest']
df_calls = study_data['calls']
df_cnv_calls = study_data['cnv_calls']
df_reviewed = study_data['reviewed_calls']
df_qc = study_data['sample_qc']

# ------- Table of cancer types -------------------------------------------------------------------
table_cancer_types(df_manifest)

# ------ QC tables --------------------------------------------------------------------------------
patients_w_fails = df_qc.loc[df_qc['total'] == 'FAIL', 'patient_id'].unique()
df_manifest_passed = df_manifest[~df_manifest['patient_id'].isin(patients_w_fails)]

# Table counting up how many samples failed QC for each test
table_qc_counts_v2(df_qc)

# Table providing more details about each QC failure
table_qc_details(df_qc)

# For safety, go ahead and exclude any calls from failed samples
df_calls = df_calls[df_calls['run_sample_id'].isin(df_manifest_passed['run_sample_id'])]
assert df_calls['run_sample_id'].nunique() == len(df_manifest_passed)

# --------- Call comparison tables ----------------------------------------------------------------
df_comps = get_call_comparisons_2test(df_manifest_passed, df_calls)
groups = compare_collection_groups(df_manifest_passed, df_comps)
combined_groups = compare_combined_collections(groups)

# Make a summary of the PPA and NPA values across all collections
# One table for PPA and one for NPA
correction_types = {1: 'raw', 2: 'raw', 3: 'raw'}
tables_ppa_npa_collections_summary(groups, correction_types)

# For comparing against acceptance criteria
ppa_352, npa_352 = tables_ppa_npa_combined_v4(groups, combined_groups, table_out_dir)

# Big table comparing the variant classes and special classes
table_call_count_summary(combined_groups, collections=[4], sections=['variant_class', 'special'],
                         correction_type='raw', ci_type='binomial', n_tests=2, add_fraction=False,
                         long_assays=True, metric_linebreak=True)

# PPA and NPA for each sample collection (variant categories)
table_call_count_summary(groups, collections=[1, 2, 3], sections=['variant_category'],
                         correction_type='raw', ci_type='binomial', n_tests=2)

# PPA and NPA for each sample collection (variant classes and special classes)
for collection in [1, 2, 3]:
    table_call_count_summary(groups, collections=[collection],
                             sections=['variant_class', 'special'],
                             correction_type='raw', ci_type='binomial', n_tests=2)

# Add a table for the performance for eaach CNV and fusion
# This information already exists in another table
# table_call_count_summary(combined_groups, collections=[4], sections=['single_gene'],
#                          correction_type='raw', ci_type='binomial', n_tests=2)

# --------- Conditioning on LDT results -----------------------------------------------------------
df_comps = get_call_comparisons(df_manifest_passed, df_calls)
groups = compare_collection_groups(df_manifest_passed, df_comps)
combined_groups = compare_combined_collections(groups)

# Show theta values and priors for variant categories
# table_prior_and_theta_values(groups)

# Show the conditioned results for variant categories for collection 3
table_call_count_summary(groups, collections=[3], sections=['variant_category'],
                         correction_type='cond', ci_type='bootstrap', n_tests=3)

# Show the conditioned results for variant classes for collection 3
table_call_count_summary(groups, collections=[3], sections=['variant_class', 'special'],
                         correction_type='cond', ci_type='bootstrap', n_tests=3)

# Describe the reviewed calls
table_variant_override_info(df_reviewed)

# --------- Focal CNVs ----------------------------------------------------------------------------
# Black sheep analysis, focal CNVs. Generate a LaTeX table for PPV and NPV
compare_focal_calls_coll123(df_calls, df_manifest_passed)

# Look at amplification status for genes on the same chromosome
compare_nonfocal_status(df_cnv_calls, patients_w_fails)

# ---------- Looking at some select discordances -------------------------------------------------
# A0132819, BRCA2, 32972742, TGTC>T
is_var = (df_comps['patient_id'] == 'A0132819') & df_comps['variant_key_nt'].str.contains('BRCA2')
df_comps[is_var]

is_var = (df_comps['patient_id'] == 'A0132834') & df_comps['variant_key_nt'].str.contains('TP53')
df_comps[is_var]

# ---------- Make a table with just the restricted reportable range -------------------------------
from D000111d.settings import PRIORS_RRR
study_data = load_study_data(
    line_data_path=os.path.join(DATA_DIR, 'line_data_new.xlsx'), use_rrr=True)
df_manifest = study_data['manifest']
df_calls = study_data['calls']
df_qc = study_data['sample_qc']

patients_w_fails = df_qc.loc[df_qc['total'] == 'FAIL', 'patient_id'].unique()
df_manifest_passed = df_manifest[~df_manifest['patient_id'].isin(patients_w_fails)]
df_calls = df_calls[df_calls['run_sample_id'].isin(df_manifest_passed['run_sample_id'])]
assert (df_calls['run_sample_id'].nunique() == 816) & (len(df_manifest_passed) == 868)

df_comps = get_call_comparisons_2test(df_manifest_passed, df_calls)

groups = compare_collection_groups(df_manifest_passed, df_comps, priors=PRIORS_RRR)
combined_groups = compare_combined_collections(groups)
ppa_353, npa_353 = tables_ppa_npa_combined_v4(groups, combined_groups, table_out_dir, suffix='_rrr')

# Make a table combining the PPA and NPA values from the two ranges for a comparsion
tables_ppa_npa_combined_bips(ppa_352, npa_352, ppa_353, npa_353)

