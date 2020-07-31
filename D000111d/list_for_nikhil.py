"""
Nikhil needs to count up how many samples that MDL has processed. I want to get our earliest
version of the manifest.
"""
import os
import pandas as pd
from D000111.settings import DATA_DIR as ORIG_DATA_DIR
from D000111b.settings import DATA_DIR as LATER_DATA_DIR


df1 = pd.read_pickle(os.path.join(ORIG_DATA_DIR, 'manifest_v1.pickle'))
df7 = pd.read_pickle(os.path.join(LATER_DATA_DIR, 'manifest_v7.pickle'))

df1_mdl = df1[df1['assay'] == 'mdl']
df7_mdl = df7[df7['assay'] == 'mdl']

df7_mdl[~df7_mdl['run_sample_id'].isin(df1_mdl['run_sample_id'])]
df1_mdl[~df1_mdl['run_sample_id'].isin(df1_mdl['run_sample_id'])]

df_needed = df7[['assay', 'run_sample_id', 'runid', 'patient_id', 'collection']]
df_wide = pd.pivot_table(df_needed, index=['collection', 'patient_id'], columns=['assay'],
                         values='run_sample_id', aggfunc=lambda x: x)
df_wide = df_wide.drop('ldt', axis=1)
df_wide = df_wide[~df_wide['mdl'].isin(['label_switch', 'missing', 'power_failure'])]
df_wide.reset_index().to_csv('~/Desktop/paired_ids.csv', index=False)

df_tall = df7.copy()
df_tall = df_tall[['run_sample_id', 'runid', 'assay', 'panel', 'biobank']].copy()
df_tall = df_tall[~df_tall['assay'].isin(['ldt'])]
df_tall = df_tall[df_tall['run_sample_id'].isin(df_wide['mdl']) |
        df_tall['run_sample_id'].isin(df_wide['ivd']) ]
df_tall.reset_index().to_csv('~/Desktop/tall_info.csv', index=False)



