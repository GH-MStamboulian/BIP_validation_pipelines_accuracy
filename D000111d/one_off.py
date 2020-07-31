"""
Get a distribution of indel lengths from production
"""
from joblib import Memory
from av360.extract.extract_calls import get_catted_call_tables
from av360.historical.identifiers import get_past_unique_ids
from settings import CACHE_DIR

memory = Memory(CACHE_DIR, verbose=0)
c_get_catted_call_tables = memory.cache(get_catted_call_tables)
c_get_past_unique_ids = memory.cache(get_past_unique_ids)


# Use av360 library to get past samples and calls
df_past = c_get_past_unique_ids(
    'ghdb_prod', start_date='2000-01-01', end_date='2019-06-24', panels=['GH2.11'],
    bip_versions=['3.5.2'], production=True, qc_pass=True, patient_info=False)
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

df_past_calls['ins_size'] = df_past_calls['ins_str'].str.len()
df_indels = df_past_calls[df_past_calls['variant_type'] == 'indel']


# Separate indels into insertions and deletions. Complex indels are double counted
df_sizes = df_indels[['runid', 'run_sample_id', 'gene', 'exon', 'ins_size', 'del_size']]
df_sizes['category'] = 'nothing'
df_sizes['size'] = -1
is_ins = df_sizes['ins_size'] > 0
is_del = df_sizes['del_size'] > 0
met_str = ((df_sizes['gene'] == 'MET') & (df_sizes['exon'] == 14)).\
    replace({True: '_met14', False: ''})

df_sizes.loc[is_ins & ~is_del, 'category'] = 'insertion' + met_str
df_sizes.loc[is_ins & ~is_del, 'size'] = df_sizes.loc[is_ins & ~is_del, 'ins_size']

df_sizes.loc[~is_ins & is_del, 'category'] = 'deletion' + met_str
df_sizes.loc[~is_ins & is_del, 'size'] = df_sizes.loc[~is_ins & is_del, 'del_size']

df_sizes.loc[is_ins & is_del, 'category'] = 'complex' + met_str
df_sizes.loc[is_ins & is_del, 'size'] = \
    df_sizes.loc[is_ins & is_del, ['ins_size', 'del_size']].max(axis=1)

# Get the min and max sizes for each category
df_sizes.groupby('category').agg({'size': [min, max, len]})
#                size
#                 min  max   len
# category
# complex           2   40   126
# complex_met14    25   25     1
# deletion          1   65  3001
# deletion_met14    1  853    27
# insertion         1   33  1010

# How many indels are greater than 30bp
df_long = df_sizes.loc[df_sizes['size'] > 30, ['runid', 'run_sample_id']].drop_duplicates()
len(df_long) / len(df_past)
sum(df_sizes['size'] > 30) / len(df_sizes)

df_sizes_small = df_sizes[~df_sizes['category'].str.contains('met14')]
df_long = df_sizes_small.loc[df_sizes_small['size'] > 30, ['runid', 'run_sample_id']].drop_duplicates()
len(df_long) / len(df_past)
sum(df_sizes_small['size'] > 30) / len(df_sizes_small)


from plotnine import ggplot, aes_string, geom_histogram, theme_minimal
# meh, save plotting for later

