"""
Place to keep all of the code for making tables in the report.
"""
import os
import re
from collections import OrderedDict
import numpy as np
import pandas as pd
from D000111d.settings import LATEX_DIR, LONG_ASSAY_NAMES, PPA_AC, NPA_AC
from av360.utils.utils import make_qc_fail_description
from av360.extract.settings import CLINICAL_CLASS_NAMES


PANEL_TYPES = {'cs': 'Clinically Relevant', 'panel': 'Panel-wide'}
VARIANT_TYPES = OrderedDict([
    ('snv', 'SNV'),
    ('indel', 'Indel'),
    ('cnv', 'CNA'),
    ('fusion', 'Fusion')
])
CLINICAL_CLASS_NAMES.\
    update({'homopolymer': 'Homopolymer', 'long_indel': 'Long Indel'})
LATEX_TABLE_DIR = os.path.join(LATEX_DIR, 'tables')


# Order the clinical class names by SNV, indel, CNA, and then fusion (cause that's what we do)
CLASS_ORDER = [
    'braf_activating_snvs', 'brca1_inactivating_snv', 'brca2_inactivating_snv',
    'egfr_l858r', 'egfr_t790m', 'kras_activating_snvs', 'nras_activating_snvs',
    'other_egfr_activating_snvs',
    'brca2_inactivating_indel', 'brca1_inactivating_indel', 'egfr_activating_indel',
    'erbb2_activating_indel', 'other_egfr_activating_indel',
    'homopolymer', 'long_indel',
    'erbb2_amplification', 'met_amplification',
    'ntrk1_fusion', 'ret_fusion', 'ros1_fusion', 'alk_fusion']
CLINICAL_CLASS_NAMES = OrderedDict([(short_name, CLINICAL_CLASS_NAMES[short_name])
                                    for short_name in CLASS_ORDER])


def ppa_formatter(ppa):
    if np.isnan(ppa):
        return 'NA'
    return '{:0.1f}'.format(100 * ppa)


def npa_formatter(npa):
    if np.isnan(npa):
        return 'NA'
    return '{:0.3f}'.format(100 * npa)


def ppa_str_formatter(ppa_str):
    return 'NA' if ('nan' in ppa_str or 'NA' in ppa_str) else ppa_str


def npa_str_formatter(npa_str):
    return 'NA' if ('nan' in npa_str or 'NA' in npa_str) else npa_str


def ci_formatter(llci, ulci, formatter, paren=True):
    out = formatter(llci) + ' - ' + formatter(ulci)
    if paren:
        out = '(' + out + ')'
    if 'nan' in out or 'NA' in out:
        out = 'NA'
    return out


def ci_stats_formatter(stats, metric, correction_type, ci_type, paren=False):
    if metric in ['ppv', 'ppa']:
        formatter = ppa_formatter
    elif metric in ['npv', 'npa', 'opa']:
        formatter = npa_formatter
    llci = stats[f'{correction_type}_{metric}_llci_{ci_type}']
    ulci = stats[f'{correction_type}_{metric}_ulci_{ci_type}']
    ci_str = ci_formatter(llci, ulci, formatter, paren=paren)
    return ci_str


def consolidate_headers(latex_str, header_rows=(2, 3), cut_column=3):
    """
    pandas annoyingly has two rows of headers when the index has names. This gets rid of those with
    some pretty hacky code.
    """
    latex_rows = latex_str.split('\n')
    top_header_row = latex_rows[header_rows[0]].split('&')
    lower_header_row = latex_rows[header_rows[1]].split('&')
    new_header_row = '&'.join(lower_header_row[0:cut_column] + top_header_row[cut_column:])
    new_rows = latex_rows[0:header_rows[0]] + [new_header_row] + latex_rows[(header_rows[1] + 1):]
    new_latex_str = '\n'.join(new_rows)
    return new_latex_str


def table_cancer_types(df_manifest):
    """
    Count up the number of cancer types in the study
    """
    df_ct = df_manifest[['patient_id', 'cancer_type']].drop_duplicates()
    df_ct_count = df_ct.groupby('cancer_type').size().reset_index().\
        rename(columns={0: 'n', 'cancer_type': 'Cancer Type'}).\
        sort_values('n', ascending=False)

    file_name = 'table_cancer_types.tex'
    latex_str = df_ct_count.to_latex(escape=False, column_format='|l|c|', index=False)
    latex_str = post_process_latex_table(latex_str, header_rows=1)
    with open(os.path.join(LATEX_TABLE_DIR, file_name), 'w') as file_pointer:
        file_pointer.write(latex_str)


def tables_ppa_npa_combined_v4(groups, combined_groups, tbl_out_dir , latex_dir=LATEX_TABLE_DIR,  suffix=''):
    """
    The third version of :py:func:`tables_ppa_npa_combined`. This does not show the
    LDT-conditioned results and it only combines collections for PPA. For NPA, only the first
    sample collection is used. It also shows the numerators and denominators since this might be
    the information that's put in the SSED.
    """
    # Just look at variant categories
    sc4_groups = {group['class']: group for group in combined_groups
                  if group['section'] == 'variant_category'}
    sc1_groups = {group['class']: group for group in groups
                  if group['section'] == 'variant_category' and group['collection'] == 1}
    assert sc4_groups.keys() == sc1_groups.keys()
    categories = list(sc4_groups.keys())

    out_tables = []
    metric_iters = zip(['ppa', 'npa'], [PPA_AC, NPA_AC], [ppa_formatter, npa_formatter],
                       [sc4_groups, sc1_groups])
    for metric, AC, formatter, sc_groups in metric_iters:
        rows = []
        for category in categories:
            panel_type, variant_type = category.split('_')
            sc_group = sc_groups[category]
            row = [('panel_type', PANEL_TYPES[panel_type]),
                   ('variant_type', VARIANT_TYPES[variant_type])]
            stats = sc_group['stats']
            cap_metric = metric.upper()
            row.extend([
                (f'Acceptable LLCI {cap_metric}', 100 * AC[(panel_type, variant_type)])
            ])
            if metric == 'ppa':
                row.extend([
                    (r'CDx+ LBP70+', stats['overall'][0]),
                    (r'CDx{\textminus} LBP70+', stats['overall'][1])
                ])
                row.extend([
                    (cap_metric, stats[f'raw_{metric}']),
                    (f'LLCI, {metric}', stats[f'raw_{metric}_llci_binomial']),
                    (f'ULCI, {metric}', stats[f'raw_{metric}_ulci_binomial']),
                ])
            else:
                row.extend([
                    (r'CDx+ LBP70{\textminus}', stats['overall'][2]),
                    (r'CDx{\textminus} LBP70{\textminus}', stats['overall'][3])
                ])
                row.extend([
                    (cap_metric, stats[f'raw_{metric}']),
                    (f'LLCI, {metric}', stats[f'raw_{metric}_llci_binomial']),
                    (f'ULCI, {metric}', stats[f'raw_{metric}_ulci_binomial']),
                ])
            row[-3:] = [(key, formatter(value)) for (key, value) in row[-3:]]
            # Bold rows if they failed AC
            # if round(stats[f'raw_{metric}_llci_bootstrap'], 1) < AC[(panel_type, variant_type)]:
            #     row[-3:] = [(key, r'\textbf{' + value + '}') for (key, value) in row[-3:]]
            rows.append(OrderedDict(row))

        df_out = pd.DataFrame(rows).\
            set_index(['panel_type', 'variant_type']).\
            sort_index()

        # Custom sort order
        product = pd.MultiIndex.from_product([
            ('Clinically Relevant', 'Panel-wide'),
            ('SNV', 'Indel', 'CNA', 'Fusion')])
        product = [p for p in product if p in df_out.index.tolist()]
        df_out = df_out.loc[product]

        df_out.index.names = [None, None]
        df_out.columns = [re.sub(', (ppa|npa)', '', col) for col in df_out.columns]
        out_tables.append(df_out)

        column_format = '|l|l|C{1.0in}|C{0.7in}|C{0.7in}|c|c|c|'
        latex_str = df_out.to_latex(escape=False, column_format=column_format, multirow=True)
        latex_str = latex_str.replace('multirow{4}{*}', 'multirow{4}{1in}')
        latex_str = post_process_latex_table(latex_str, header_rows=1)
        df_out.to_csv(tbl_out_dir+ f'{metric}_combined4{suffix}.tsv', sep = '\t')
        with open(os.path.join(latex_dir, f'{metric}_combined4{suffix}.tex'), 'w') as file_pointer:
            file_pointer.write(latex_str)

    return out_tables


def table_qc_counts(df_qc):
    """
    Make a table counting up how many samples failed or passed QC for each test.

    Nothing is returned as output. Instead, a LaTeX file is written in :py:data:`LATEX_TABLE_DIR`

    Args:
        df_qc (pandas.DataFrame): the QD data file
    """
    # Pivot tables to count the number of failing and passing samples
    df_paired = df_qc.pivot_table(
        index=['patient_id', 'collection'], columns='assay', values='total', aggfunc=lambda x: x).\
        drop('ldt', axis=1).reset_index()
    df_counts = df_paired.pivot_table(index='collection', columns=['mdl', 'cdx'],
                                      values='patient_id', aggfunc=len, fill_value=0,
                                      margins=True, margins_name='Total')

    # Formatting for LaTeX
    df_counts.columns.names = [r'\multicolumn{1}{|r|}{LPB70:}', r'\multicolumn{1}{|r|}{G360 CDx:}']
    df_counts.index.name = 'SC'
    df_counts.columns = df_counts.columns.set_levels(['{}', 'PASS', 'FAIL'], level=1)
    latex_str = df_counts.to_latex(
        multicolumn=True, multicolumn_format='|c|', escape=False, column_format='|c|c|c|c|')
    latex_str = post_process_latex_table(latex_str, header_rows=3)
    with open(os.path.join(LATEX_TABLE_DIR, 'qc_table_count.tex'), 'w') as file_pointer:
        file_pointer.write(latex_str)


def table_qc_counts_v2(df_qc):
    """
    Make a table counting up how many samples failed or passed QC for each test.

    Nothing is returned as output. Instead, a LaTeX file is written in :py:data:`LATEX_TABLE_DIR`

    Args:
        df_qc (pandas.DataFrame): the QD data file
    """
    # Pivot tables to count the number of failing and passing samples
    df_paired = df_qc.pivot_table(
        index=['patient_id', 'collection'], columns='assay', values='total', aggfunc=lambda x: x).\
        drop('ldt', axis=1).reset_index()
    df_counts = df_paired.pivot_table(index='collection', columns=['mdl', 'cdx'],
                                      values='patient_id', aggfunc=len, fill_value=0,
                                      margins=True, margins_name='Total')
    df_counts = df_counts.reset_index()
    df_counts['collection'] = ['1', '2', '3', 'Total Samples']
    df_counts.columns = ['SC',
                         'LBP70 QC Fail, G360 CDx QC Pass',
                         'LBP70 QC Pass, G360 CDx QC Pass',
                         'Total Samples']

    # Formatting for LaTeX
    latex_str = df_counts.to_latex(
        multicolumn=True, multicolumn_format='|c|', escape=False,
        column_format='|l|C{1.5in}|C{1.5in}|C{1in}|', index=False)
    latex_str = post_process_latex_table(latex_str, header_rows=1)
    with open(os.path.join(LATEX_TABLE_DIR, 'qc_table_count.tex'), 'w') as file_pointer:
        file_pointer.write(latex_str)


def find_nth(substr, string, n):
    """ Find the nth substr in string. No error checking. """
    start_loc = -1
    for _ in range(n):
        start_loc = string.find(substr, start_loc + 1)
    return start_loc


def post_process_latex_table(latex_str, header_rows=0):
    """
    Take the output of pandas ``to_latex`` function and format it to a version that converts more
    easily to a MS Word table.
    """
    # latex_str = df_counts.to_latex(
    #     multicolumn=True, multicolumn_format='|c|', escape=False, column_format='|c|c|c|c|')
    # Get rid of the booktabs weirdness
    latex_str = re.sub('toprule', r'hline', latex_str)
    latex_str = re.sub(r'\\bottomrule', '', latex_str)
    latex_str = re.sub(r'\\midrule', '', latex_str)
    # Remove any empty rows
    latex_str = '\n'.join([row for row in latex_str.split('\n') if row != ''])
    # Remove any rows with cline (a hack)
    latex_str = '\n'.join([row for row in latex_str.split('\n') if 'cline' not in row])

    # Add horizontal lines for every row (or clines if there is a multirow)
    lines = latex_str.split('\n')
    n_columns = lines[header_rows + 2].count('&') + 1
    multirows_left = [0] * n_columns
    for i_line in range(2, len(lines)-1):
        line = lines[i_line]
        if 'multirow' in line:
            cells = line.split('&')
            for i_cell, cell in enumerate(cells):
                finds = re.findall('multirow{([0-9])}', cell)
                assert len(finds) <= 1
                if finds:
                    multirows_left[i_cell] = int(finds[0])
            # Only left cells can be merged
            assert sorted(multirows_left, reverse=True) == multirows_left
            left_cline = 1 + len([x for x in multirows_left if x > 1])
            cline_str = str(left_cline) + '-' + str(n_columns)
            lines[i_line] = re.sub(r'\\\\', r'\\\\ \\cline{' + cline_str + '}', line)
            multirows_left = [x-1 for x in multirows_left]
        elif any([x > 1 for x in multirows_left]):
            left_cline = 1 + len([x for x in multirows_left if x > 1])
            cline_str = str(left_cline) + '-' + str(n_columns)
            lines[i_line] = re.sub(r'\\\\', r'\\\\ \\cline{' + cline_str + '}', line)
            multirows_left = [x-1 for x in multirows_left]
        else:
            last_position = line.rfind(r'\\')
            lines[i_line] = line[:last_position] + r'\\ \hline' + line[(last_position + 2):]
    latex_str = '\n'.join(lines)

    # Make the header rows bold text with a gray background
    lines = latex_str.split('\n')
    for i_row in range(2, 2 + header_rows):
        # We assume that there are two linese of text before the tabular content
        header_line = lines[i_row]
        line_cells = header_line.split('&')
        # All but the last column are simply surrounded by \textbf{}
        for i_col in range(len(line_cells) - 1):
            cell = line_cells[i_col]
            if cell.strip() == '' or cell.strip() == '{}':
                continue
            elif 'multicolumn' in cell:
                # Hacky. Find the third opening and closing braces
                left = find_nth('{', cell, 3)
                right = find_nth('}', cell, 3)
                middle_str = r'\textbf{' + cell[(left+1):right] + '}'
                line_cells[i_col] = cell[0:left+1] + middle_str + cell[right:]
            else:
                line_cells[i_col] = r'\textbf{' + cell.strip() + '}'
        # Make sure not to enclose the backslahses and hlines in textbf for the last column
        last_cell = line_cells[-1]
        if 'multicolumn' in last_cell:
            left = find_nth('{', last_cell, 3)
            right = find_nth('}', last_cell, 3)
            middle_str = r'\textbf{' + last_cell[(left+1):right] + '}'
            line_cells[-1] = last_cell[0:left+1] + middle_str + last_cell[right:]
        else:
            line_cells[-1] = re.sub(r'\\\\', r'}\\\\', line_cells[-1])
            line_cells[-1] = r'\textbf{' + line_cells[-1]

        header_line = ' & '.join(line_cells)
        # Make the header line have a gray background
        lines[i_row] = r'\rowcolor[gray]{.85}' + header_line
    latex_str = '\n'.join(lines)
    # print(latex_str)
    return latex_str


def table_qc_details(df_qc):
    """
    Make a table providing a string to describe each QC failure.

    Nothing is returned as output. Instead, a LaTeX file is written in :py:data:`LATEX_TABLE_DIR`

    Args:
        df_qc (pandas.DataFrame): the QD data file
    """
    status_cols = [col for col in df_qc.columns if col.endswith('_status')]
    value_cols = [re.sub('_status$', '', col) for col in status_cols]
    values_cols = [col for col in value_cols if col in df_qc.columns]
    id_cols = [col for col in df_qc.columns if col not in status_cols + value_cols]
    df_status = df_qc[id_cols + status_cols]
    df_status.columns = [re.sub('_status$', '', col) for col in df_status.columns]
    df_values = df_qc[id_cols + value_cols]

    exclude_metrics = ('total', 'sample_autoqc_total', 'sample_germline_contamination')
    qc_strings = []
    for i_row, row_status in df_status.iterrows():
        row_values = df_values.loc[i_row].drop(['sample_autoqc_total', 'total'])
        qc_strings.append(make_qc_fail_description(row_values, row_status, exclude_metrics))
    df_status['qc_description'] = qc_strings

    # Make sure that everything that fails QC has a description
    assert not any((df_status['qc_description'] == '') & (df_status['total'] == 'FAIL'))

    # Formatting for LaTeX
    df_status['qc_description'] = df_status['qc_description'].\
        str.replace('exclude_reason: ', '').\
        str.replace('correct_run_sample_id', 'Incorrect Sample ID')
    df_latex = df_status.loc[df_status['total'] == 'FAIL',
                             ['collection', 'assay', 'run_sample_id', 'qc_description']].copy()
    df_latex['assay'] = df_latex['assay'].replace(LONG_ASSAY_NAMES)
    df_latex['run_sample_id'] = df_latex['run_sample_id'].\
        replace({'qc_fail': '', 'power_failure': ''})
    df_latex = df_latex.sort_values(['collection', 'assay', 'qc_description'])
    df_latex = df_latex.rename(
        columns={'collection': 'SC', 'assay': 'Test', 'run_sample_id': 'Sample ID',
                 'qc_description': 'QC Failures'})
    latex_str = df_latex.to_latex(index=False, column_format='|r|l|l|l|')
    latex_str = post_process_latex_table(latex_str, header_rows=1)
    with open(os.path.join(LATEX_TABLE_DIR, 'qc_reason_table.tex'), 'w') as file_pointer:
        file_pointer.write(latex_str)


def table_call_count_summary(groups, collections, sections,
                             correction_type='raw', ci_type='binomial',
                             n_tests=None, add_fraction=False, long_assays=False,
                             metric_linebreak=False, latex_dir=LATEX_TABLE_DIR):
    """
    Summarize multiple contingency tables into a single table

    Args:
        groups (list): the objects containing all the information for each comparison
        collection (list): which collections should be compared (example: [1, 2, 3])
        sections (list): which "sections" should be compared (classifier for groups)
        correction_type (str): "raw" or "cond" (for conditioned)
        ci_type (str): either "binomial" or "bootstrap"
        n_tests (int): either 2 for 2-test comparison or 3 for 3-test comparison
        add_fraction (bool): if ``True``, then also show the fraction before PPA and NPA values
        long_assays (bool): if ``True``, then don't just show the letters a, b, c, and d (etc.)
            that indicate the location in the contingency tables. Instead, make the column header
            more descriptive.
    """
    coll_groups = [group for group in groups
                   if group['collection'] in collections and
                   group['section'] in sections]
    rows = []
    for group in coll_groups:
        stats = group['stats']
        if group['section'] == 'variant_category':
            panel_type, variant_type = group['class'].split('_')
            row = [('SC', group['collection']),
                   ('panel_type', PANEL_TYPES[panel_type]),
                   ('variant_type', VARIANT_TYPES[variant_type])]
            index_columns = ['SC', 'panel_type', 'variant_type']
        elif group['section'] in ['variant_class', 'special']:
            row = [('SC', group['collection']),
                   ('class', CLINICAL_CLASS_NAMES[group['class']])]
            index_columns = ['SC', 'class']
        elif group['section'] == 'single_gene':
            variant_type, gene = group['class'].split(':')
            row = [('SC', group['collection']),
                   ('variant_type', VARIANT_TYPES[variant_type]),
                   ('gene', gene)]
            index_columns = ['SC', 'variant_type', 'gene']
        else:
            assert False
        if len(collections) == 1:
            row = row[1:]
            index_columns = index_columns[1:]

        # How many tests are there?
        if n_tests is None:
            n_tests = group['n_tests']
        overall = stats['overall']

        # Strings for PPA and NPA
        metric_strings = []
        for metric, formatter in zip(['ppa', 'npa'], [ppa_formatter, npa_formatter]):
            metric_ci = ci_type
            if isinstance(ci_type, dict):
                metric_ci = ci_type[metric]
            if metric_linebreak:
                metric_str = \
                    r'\makecell{' + \
                    formatter(stats[f'{correction_type}_{metric}']) + \
                    r' \\ (' + \
                    formatter(stats[f'{correction_type}_{metric}_llci_{metric_ci}']) + \
                    ' - ' + \
                    formatter(stats[f'{correction_type}_{metric}_ulci_{metric_ci}']) + \
                    ')}'
            else:
                metric_str = \
                    formatter(stats[f'{correction_type}_{metric}']) + \
                    ' (' + \
                    formatter(stats[f'{correction_type}_{metric}_llci_{metric_ci}']) + \
                    ' - ' + \
                    formatter(stats[f'{correction_type}_{metric}_ulci_{metric_ci}']) + \
                    ')'

            if add_fraction:
                assert correction_type == 'raw'
                assert ci_type == 'binomial'
                assert n_tests == 2
                if metric == 'ppa':
                    num, denom = overall[0], overall[0] + overall[1]
                else:
                    num, denom = overall[3], overall[2] + overall[3]
                metric_str = r'$\sfrac{' + str(num) + '}{' + str(denom) + '} = $' + metric_str
            metric_strings.append(metric_str)
        ppa_str, npa_str = metric_strings

        # Show the number of concordances and discordances
        column_formats = ['l'] * len(row)
        if n_tests == 2:
            if len(overall) > 4:
                a0, b0, c0, d0, a1, b1, c1, d1 = overall
                overall = [a0 + a1, b0 + b1, c0 + c1, d0 + d1]
            row.extend([
                ('CDx+ LBP70+' if long_assays else '$a$', overall[0]),
                (r'CDx{\textminus} LBP70+' if long_assays else '$b$', overall[1]),
                (r'CDx+ LBP70{\textminus}' if long_assays else '$c$', overall[2]),
                (r'CDx{\textminus} LBP70{\textminus}' if long_assays else '$d$', overall[3]),
            ])
            count_column_formats = ['l'] * 4
            if long_assays:
                count_column_formats = ['C{0.54in}', 'C{0.54in}', 'C{0.54in}', 'C{0.54in}']
        elif n_tests == 3:
            row.extend([
                ('$a_0$', overall[0]),
                ('$b_0$', overall[1]),
                ('$c_0$', overall[2]),
                ('$d_0$', overall[3]),
                ('$a_1$', overall[4]),
                ('$b_1$', overall[5]),
                ('$c_1$', overall[6]),
                ('$d_1$', overall[7]),
            ])
            count_column_formats = ['l'] * 8
        elif n_tests == 0:
            count_column_formats = []
        row.extend([
            ('PPA', ppa_str_formatter(ppa_str)),
            ('NPA', npa_str_formatter(npa_str))
        ])
        if not metric_linebreak:
            column_formats += count_column_formats + ['l', 'l']
        else:
            column_formats += count_column_formats + ['c', 'c']

        assert len(column_formats) == len(row)
        rows.append(OrderedDict(row))

    df_out = pd.DataFrame(rows).set_index(index_columns)
    if 'variant_category' in sections:
        # Sort in order of snv, indel, cna, then fusion
        if len(collections) > 1:
            product = pd.MultiIndex.from_product([
                (1, 2, 3), ('Clinically Relevant', 'Panel-wide'),
                ('SNV', 'Indel', 'CNA', 'Fusion')])
        else:
            product = pd.MultiIndex.from_product([
                ('Clinically Relevant', 'Panel-wide'),
                ('SNV', 'Indel', 'CNA', 'Fusion')])
        product = [p for p in product if p in df_out.index.tolist()]
        df_out = df_out.loc[product]
        df_out.index.names = ['{}' if x != 'SC' else x for x in index_columns]
        # Make the TeX
        column_formats[1] = 'L{0.7in}'
        column_format = '|' + '|'.join(column_formats) + '|'
        latex_str = df_out.to_latex(escape=False, column_format=column_format, multirow=True)
        # A hardcoded hack for the "Panel Type" column
        latex_str = latex_str.replace(r'\multirow{4}{*}', r'\multirow{4}{0.7in}')
        latex_str = latex_str.replace(r'\multirow{2}{*}', r'\multirow{2}{0.7in}')
        new_latex_str = consolidate_headers(latex_str, cut_column=len(index_columns))
    elif 'variant_class' in sections or 'special' in sections:
        # Sort in the order defined in the CLINICAL_CLASS_NAMES variable
        df_out = df_out.loc[CLINICAL_CLASS_NAMES.values()]
        df_out.index.names = ['{}' if x != 'SC' else x for x in index_columns]
        # Make the TeX
        column_format = '|' + '|'.join(column_formats) + '|'
        latex_str = df_out.to_latex(escape=False, column_format=column_format, multirow=True)
        new_latex_str = consolidate_headers(latex_str, cut_column=len(index_columns))
    elif 'single_gene' in sections:
        row_order = [('CNA', 'MET'), ('CNA', 'ERBB2'), ('Fusion', 'ALK'), ('Fusion', 'RET'),
                     ('Fusion', 'ROS1'), ('Fusion', 'NTRK1')]
        df_out = df_out.loc[row_order]
        column_format = '|' + '|'.join(column_formats) + '|'
        latex_str = df_out.to_latex(escape=False, column_format=column_format, multirow=True)
        new_latex_str = consolidate_headers(latex_str, cut_column=len(index_columns))
    else:
        assert False

    collections_str = ''.join([str(c) for c in collections])
    sections_str = ''.join(sections)
    if isinstance(ci_type, dict):
        ci_str = '_'.join([str(metric) + ':' + str(metric_ci)
                           for (metric, metric_ci) in ci_type.items()])
    else:
        ci_str = ci_type
    file_name = f'coll{collections_str}_{sections_str}_{correction_type}_{ci_str}'
    file_name = f'{file_name}_{n_tests}test.tex'
    new_latex_str = post_process_latex_table(new_latex_str, header_rows=1)
    with open(os.path.join(latex_dir, file_name), 'w') as file_pointer:
        file_pointer.write(new_latex_str)


def tables_ppa_npa_collections_summary(groups, correction_types, latex_dir=LATEX_TABLE_DIR):
    """
    Make a summary of the PPA and NPA values across all collections. One table for PPA and one for
    NPA

    Args:
        correction_types (dict): key is the collection and value is the type of correction that was
            done for PPA or NPA calculation
    """
    coll_groups = [group for group in groups
                   if group['section'] == 'variant_category']

    for metric, formatter in zip(['ppa', 'npa'], [ppa_formatter, npa_formatter]):
        rows = []
        # Make a tall table and then widen it
        for group in coll_groups:
            stats = group['stats']
            ref_calls = stats['ref_pos'] if metric == 'ppa' else stats['ref_neg']
            panel_type, variant_type = group['class'].split('_')
            correction_type = correction_types[group['collection']]
            rows.append({
                'variant_type': VARIANT_TYPES[variant_type],
                'panel_type': PANEL_TYPES[panel_type],
                metric.upper(): formatter(stats[f'{correction_type}_{metric}']),
                'n': ref_calls,
                'collection': 'Collection ' + str(group['collection'])
            })
        df_tall = pd.DataFrame(rows)
        df_taller = pd.melt(df_tall, id_vars=['collection', 'panel_type', 'variant_type'],
                            var_name='metric')
        df_out = df_taller.\
            pivot_table(index=['panel_type', 'variant_type'],
                        columns=['collection', 'metric'],
                        values='value', aggfunc=lambda x: x)

        # Custom sort order
        product = pd.MultiIndex.from_product([
            ('Clinically Relevant', 'Panel-wide'),
            ('SNV', 'Indel', 'CNA', 'Fusion')])
        product = [p for p in product if p in df_out.index.tolist()]
        df_out = df_out.loc[product]

        df_out.index.names = [None, None]
        df_out.columns.names = [None, None]
        cor_types = ''.join(pd.Series(correction_types).sort_index())
        file_name = f'{metric}_summary_collections_{cor_types}.tex'
        latex_str = df_out.to_latex(
            escape=False, column_format='|l|l|c|c|c|c|c|c|', multirow=True,
            multicolumn_format='c|')
        latex_str = post_process_latex_table(latex_str, header_rows=2)
        with open(os.path.join(latex_dir, file_name), 'w') as file_pointer:
            file_pointer.write(latex_str)


def table_prior_and_theta_values(groups):
    """ Display the priors and theta values """
    theta_groups = [group for group in groups
                    if group['n_tests'] == 3 and group['section'] == 'variant_category']
    rows = []
    for group in theta_groups:
        stats = group['stats']
        panel_type, variant_type = group['class'].split('_')
        row = OrderedDict([
            ('panel_type', PANEL_TYPES[panel_type]),
            ('variant_type', VARIANT_TYPES[variant_type]),
            ('Unique Sites', group['n_vars']),
            ('Average Variants Per Sample', '{:0.3f}'.format(group['var_per_sample'])),
            (r'$P(\text{LDT}^+)$', '{:0.6f}'.format(group['prior_freq'])),
            (r'$\theta$', '{:0.2f}'.format(stats['theta']))
        ])
        rows.append(row)
    df_out = pd.DataFrame(rows).\
        set_index(['panel_type', 'variant_type']).\
        sort_index()
    df_out.index.names = [None, None]
    file_name = 'prior_theta_values_variant_categories.tex'
    with open(os.path.join(LATEX_TABLE_DIR, file_name), 'w') as file_pointer:
        df_out.to_latex(file_pointer, escape=False,
                        column_format='llC{0.75in}C{1in}C{0.75in}C{0.5in}')


def table_variant_override_info(df_reviewed):
    """
    Make a table explaining all the reasons for why some variant calls were changed.

    Args:
        df_reviewed (pandas.DataFrame): output of ``get_reviewed_calls()``
    """
    short_fusion_suppress_note = 'Near-fusion suppression'
    short_snv_alignment_note = 'Indel realignment, suppressed SNV'
    short_indel_alignment_note = 'Different indel alignments'
    short_complex_indel_note = "Equivalent complex indels"
    notes = {('A0132572', 'TP53.17.7579371.TGCC>AGCCT'): short_fusion_suppress_note,
             ('A0132572', 'TP53.17.7579476.G>T'): short_fusion_suppress_note,
             ('A0132581', 'NOTCH1.9.139412214.G>A'): short_fusion_suppress_note,
             ('A0133866', 'ERBB2.17.37881359.GT>G'): short_fusion_suppress_note,
             ('A0132826', 'ERBB2.17.37880996.T>A'): short_snv_alignment_note,
             ('A0132827', 'ERBB2.17.37880996.T>A'): short_snv_alignment_note,
             ('A0132828', 'ERBB2.17.37880996.T>A'): short_snv_alignment_note,
             ('A0133856', 'ERBB2.17.37880998.G>T'): short_snv_alignment_note,
             ('A0132834', 'TP53.17.7579702.G>A'): short_snv_alignment_note,
             ('A0132815', 'BRCA2.13.32914477.CGCAAGACAAGT>C'): short_complex_indel_note,
             ('A0132815', 'BRCA2.13.32914478.GCAAGACAAGTGTT>ATC'): short_complex_indel_note,
             ('A0132817', 'BRCA2.13.32914314.AAGTTTCTAAAATATCACCTTGTGAT>TATCA'):
             short_complex_indel_note,
             ('A0132817', 'BRCA2.13.32914314.AAGTTTCTAAAATATCACCTTGTGAT>TATTA'):
             short_complex_indel_note,
             ('A0132816', 'BRCA2.13.32929220.TAAAACTAAA>T'): short_complex_indel_note,
             ('A0132816', 'BRCA2.13.32929220.TAAAACTAAAT>GTTC'): short_complex_indel_note,
             ('A0132872', 'APC.5.112174757.GA>G'): short_indel_alignment_note,
             ('A0132834', 'TP53.17.7579685.AACCCTTGTCCTTACCAG>A'): short_indel_alignment_note,
             ('A0132819', 'BRCA2.13.32972742.TGTC>T'): short_indel_alignment_note}

    pa_call_to_call = {0: '-', 1: '+', 2: '+', 'remove': '-',
                       '0': '-', '1': '+', '2': '+'}
    rows = []
    for _, row in df_reviewed.iterrows():
        var_key = row['variant_key_nt']
        patient_id = row['patient_id']
        gene, chrom, position, mut_nt = var_key.split('.')
        mdl_original_call = pa_call_to_call[row['mdl']]
        cdx_original_call = pa_call_to_call[row['cdx']]
        mdl_new_call = pa_call_to_call[row['mdl_reviewed']]
        cdx_new_call = pa_call_to_call[row['cdx_reviewed']]
        new_row = OrderedDict([
            ('Patient', patient_id),
            ('Type', row['variant_type']),
            ('Gene', gene),
            ('Position', int(position)),
            ('Mutation', mut_nt),
            ('LBP70', f'{mdl_original_call} ({mdl_new_call})'),
            ('CDx', f'{cdx_original_call} ({cdx_new_call})'),
            ('Reason', notes.get((patient_id, var_key), ''))
        ])
        rows.append(new_row)

    df_out = pd.DataFrame(rows)
    df_out = df_out.sort_values(['Reason', 'Patient'])
    df_out['Type'] = df_out['Type'].replace(VARIANT_TYPES)
    df_out = df_out.set_index(['Reason', 'Patient', 'Type', 'Gene', 'Position'])
    column_format = '|L{1.1in}|l|l|l|l|L{1.6in}|l|l|'
    latex_str = df_out.to_latex(index=True, escape=True, multirow=True,
                                column_format=column_format)
    new_latex_str = consolidate_headers(latex_str, cut_column=df_out.index.nlevels)
    for reason in set(notes.values()):
        new_latex_str = new_latex_str.replace('{*}{' + reason, '{1.1in}{' + reason)
    file_name = 'reviewed_calls.tex'
    new_latex_str = post_process_latex_table(new_latex_str, header_rows=1)
    with open(os.path.join(LATEX_TABLE_DIR, file_name), 'w') as file_pointer:
        file_pointer.write(new_latex_str)


def tables_ppa_npa_combined_bips(ppa_352, npa_352, ppa_353, npa_353):
    """
    Combine the PPA and NPA tables from the different BIP versions. The inputs are data frames that
    come from :py:func:`tables_ppa_npa_combined_v4`
    """
    sum_tables = {'PPA': [ppa_352, ppa_353], 'NPA': [npa_352, npa_353]}
    for metric, metric_tables in sum_tables.items():
        table_352, table_353 = metric_tables
        df_out = pd.concat([table_352, table_353], axis=0, keys=['3.5.2', '3.5.3']).\
            reorder_levels([1, 2, 0])

        # Custom sort order
        product = pd.MultiIndex.from_product([
            ('Clinically Relevant', 'Panel-wide'),
            ('SNV', 'Indel', 'CNA', 'Fusion'),
            ('3.5.2', '3.5.3')]).tolist()
        product = [p for p in product if p in df_out.index.tolist()]
        df_out = df_out.loc[product]

        # Change the indeces
        df_out = df_out.reset_index().\
            rename(columns={'level_0': 'CR', 'level_1': 'Variant Category', 'level_2': 'BIP'}).\
            set_index(['CR', 'Variant Category', f'Acceptable LLCI {metric}', 'BIP'])
        df_out.index.names = ['', '', f'Acceptable LLCI {metric}', 'BIP']

        # Format and make latex string
        column_format = '|l|l|L{1.0in}|L{0.5in}|C{0.7in}|C{0.7in}|c|c|c|'
        latex_str = df_out.to_latex(escape=False, column_format=column_format, multirow=True)
        latex_str = latex_str.replace('multirow{4}{*}', 'multirow{4}{0.75in}')
        latex_str = latex_str.replace('multirow{8}{*}', 'multirow{8}{0.75in}')
        latex_str = consolidate_headers(latex_str, header_rows=(2, 3), cut_column=4)
        latex_str = post_process_latex_table(latex_str, header_rows=1)
        file_path = os.path.join(LATEX_TABLE_DIR, f'{metric}_combined_bips.tex')
        with open(file_path, 'w') as file_pointer:
            file_pointer.write(latex_str)
