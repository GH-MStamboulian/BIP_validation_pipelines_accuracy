""" Custom analyses for focal CNVs """
import os
import pandas as pd
from statsmodels.stats.proportion import proportion_confint
from D000111d.tables import consolidate_headers, ppa_formatter, npa_formatter, ci_formatter, \
    post_process_latex_table
from D000111d.settings import LATEX_DIR, LONG_ASSAY_NAMES


LATEX_TABLE_DIR = os.path.join(LATEX_DIR, 'tables')


def compare_focal_calls_coll123(df_calls, df_manifest_passed):
    """
    Don't separate the collections. Add them all together.
    """
    # Get just the information that's needed
    df_cnv = df_calls.loc[df_calls['variant_type'] == 'cnv',
                          ['patient_id', 'assay', 'variant_key_nt', 'call', 'collection']].copy()
    df_cnv['call'] = df_cnv['call'].astype(int)
    df_cnv_comp = pd.pivot_table(df_cnv, index=['collection', 'patient_id', 'variant_key_nt'],
                                 columns='assay', values='call',
                                 fill_value=0, aggfunc=lambda x: x).reset_index()

    # Convert calls to 0=not detected, 1=detected, 2=reported
    df_cnv_comp['cdx'] = df_cnv_comp['cdx'].replace({0: 0, 1: 1, 2: 2, 3: 1})
    df_cnv_comp['ldt'] = df_cnv_comp['ldt'].replace({1: 2})
    df_cnv_comp['mdl'] = df_cnv_comp['mdl'].replace({1: 2})

    # 2-test comparisons
    n_unique_vars = 2
    n_patients = df_manifest_passed['patient_id'].nunique()
    n_vars = n_patients * n_unique_vars
    df2 = pd.pivot_table(df_cnv_comp[df_cnv_comp['collection'].isin([1, 2, 3])],
                         index=['cdx'], columns='mdl', aggfunc=len,
                         values='patient_id', fill_value=0)
    n_pos = df2.sum().sum() - df2.loc[0, 0]
    n_neg = n_vars - n_pos
    df2.loc[0, 0] = n_neg

    # Calculate PPV and NPV
    true_pos_calls = df2.loc[2, 2]
    pos_calls = df2.loc[2, :].sum()
    true_neg_calls = df2.loc[0:1, 0].sum()
    neg_calls = df2.loc[0:1, :].sum().sum()
    ppv = true_pos_calls / pos_calls
    npv = true_neg_calls / neg_calls

    # And the confidence intervals
    cis = proportion_confint(true_pos_calls, pos_calls, method='beta')
    ci_str = ci_formatter(cis[0], cis[1], ppa_formatter)
    ppv_str = ppa_formatter(ppv) + ' ' + ci_str
    cis = proportion_confint(true_neg_calls, neg_calls, method='beta')
    ci_str = ci_formatter(cis[0], cis[1], npa_formatter)
    npv_str = npa_formatter(npv) + ' ' + ci_str

    # Insert PPV and NPV into the table
    df2['metric'] = ''
    df2.loc[0, 'metric'] = 'NPV = ' + npv_str + r' \% '
    df2.loc[2, 'metric'] = 'PPV = ' + ppv_str + r' \% '

    # Do some formatting
    df2 = df2.rename(index={0: 'ND, CDx', 1: 'D, CDx', 2: 'R, CDx'}, level='cdx')
    df2 = df2.rename(columns={0: r'LBP70\textminus', 2: 'LBP70+', 'metric': '{}'})
    df2.columns.name = None
    df2.index.names = [None]
    latex_str = df2.to_latex(escape=False)
    latex_str = post_process_latex_table(latex_str, header_rows=2)

    # and write to a LaTeX table
    with open(os.path.join(LATEX_TABLE_DIR, 'focal_cnvs_coll123.tex'), 'w') as file_pointer:
        file_pointer.write(latex_str)


def compare_nonfocal_status(df_cnv_calls, patients_w_fails):
    """
    For CDx, nonfocal CNV calls, give the calls for all othere genes on the chromosome (both the
    MDL and the CDx calls)
    """
    df_cnv = df_cnv_calls[~df_cnv_calls['patient_id'].isin(patients_w_fails)]
    neighbors = {
        # On 17q (TP53 and CHD3 on 17p)
        'ERBB2': ['NF1', 'CDK12', 'BRCA1', 'XYLT2', 'TP53', 'CHD3'],
        # On 7q (EGFR on 7p)
        'MET': ['CDK6', 'IMPDH1', 'SMO', 'BRAF', 'EZH2', 'RHEB', 'EGFR']
    }
    for gene, neighbor_genes in neighbors.items():
        nf_patients = df_cnv.loc[(df_cnv['call'] == 3) &
                                 (df_cnv['variant_key_nt'] == gene), 'patient_id']
        df_neighbors = df_cnv[df_cnv['patient_id'].isin(nf_patients) &
                              df_cnv['variant_key_nt'].isin(neighbor_genes + [gene]) &
                              (df_cnv['assay'] != 'ldt')]
        pdf = pd.pivot_table(df_neighbors, index='patient_id', columns=['variant_key_nt', 'assay'],
                             values='call', aggfunc=lambda x: x)
        n_genes = pdf.columns.get_level_values('variant_key_nt').nunique()
        pdf.loc[:, (slice(None), 'mdl')] = \
            pdf.loc[:, (slice(None), 'mdl')].replace({0: r'\textminus', 1: '+'})
        pdf.loc[:, (slice(None), 'cdx')] = \
            pdf.loc[:, (slice(None), 'cdx')].replace({0: r'\textminus', 1: 'D', 2: 'F', 3: 'D'})
        pdf = pdf.rename(columns=LONG_ASSAY_NAMES).reset_index().\
            rename(columns={'patient_id': 'Patient'})
        column_format = '|l' + '|c|c' * n_genes + '|'
        latex_str = pdf.to_latex(index=False, multicolumn_format='c|', multicolumn=True,
                                 escape=False, column_format=column_format)
        latex_str = post_process_latex_table(latex_str, header_rows=2)
        with open(os.path.join(LATEX_TABLE_DIR, f'nonfocal_cnvs_{gene}.tex'), 'w') as file_pointer:
            file_pointer.write(latex_str)
