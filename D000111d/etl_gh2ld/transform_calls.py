"""
The number of mutation and filtration steps for the call data was complex enough that it deserved
it's own file.
"""
import os
import pandas as pd
import pysam
from av360.transform.transform_calls import is_av_variant
from av360.backend.find_backend import start_backend
from pybedtools import BedTool
from D000111d.settings import BED_FILE_PATHS, BLACKLIST_BED_FILE_PATHS


HG19_PATH = '{HOME}/.cache/genomes/hg19.fa'.format(HOME=os.environ['HOME'])
hg19 = pysam.FastaFile(HG19_PATH)


def transform_calls(df_calls, df_reviewed, df_manifest):
    """
    Wrapper for all of the smaller transformation functions
    """
    # Only look at calls that are in the manifest
    df_calls = df_calls[df_calls['run_sample_id'].isin(df_manifest['run_sample_id'])]

    # Only look at CNVs and fusions that we are validating for CDx
    df_calls = df_calls[is_av_variant(df_calls)]

    # Add some special variant classifications relevant to just this study
    df_calls['homopolymer'] = annotate_homopolymers(df_calls)
    df_calls['long_indel'] = annotate_long_indels(df_calls)

    # Update the manually adjudicated calls
    df_calls = override_calls(df_calls, df_reviewed, df_manifest)

    # Remove any fusion and CNV calls where call = 0
    df_calls = df_calls[(~df_calls['variant_type'].isin(['fusion', 'cnv']) |
                         (df_calls['call'] != 0))]

    # Add columns that annotate if the call is in any of the panels in the manifest
    df_calls = pd.merge(df_calls, df_manifest[['run_sample_id', 'panel', 'assay']], 
                        how='left')
    df_reportable = get_reportable_status(df_calls, 'variant_key_nt')
    df_reportable.columns = ['panel_' + col for col in df_reportable.columns]
    panel_columns = df_reportable.columns
    df_reportable = df_reportable.reset_index()
    df_calls = pd.merge(df_calls, df_reportable, how='left')

    # Apply the indel filters described in the protocol
    var_comment = df_calls['variant_comment'].str
    exclude_indel_filter = (df_calls['variant_type'] == 'indel') & \
        (var_comment.contains('oncogene fs var_ds <= 1').fillna(False) |
         var_comment.contains('single molecule indel present').fillna(False))
    df_calls = df_calls[~exclude_indel_filter]

    # For MDL samples, we need to override the "pa_call" to be 2 if call > 0
    is_mdl_cnv = (df_calls['assay'] == 'mdl') & (df_calls['variant_type'] == 'cnv')
    df_calls.loc[is_mdl_cnv, 'pa_call'] = \
        df_calls.loc[is_mdl_cnv, 'call'].replace({1: 2})

    # Minimize the number of columns that are returned
    small_columns = \
        ['runid', 'run_sample_id', 'call', 'pa_call', 'variant_type', 'clinical_class',
         'variant_key_nt', 'percentage', 'mut_fams', 'fams', 'mut_nt', 'mut_aa', 'position',
         'gene', 'homopolymer', 'long_indel', 'splice_effect', 'chrom', 'del_size', 'in_rrr',
         'rm_reportable'] + \
        panel_columns.tolist()
    small_columns = [col for col in small_columns if col in df_calls]
    df_small_calls = df_calls[small_columns]

    return df_small_calls


def transform_cnv_calls(df_calls):
    """
    The focal CNV analysis requires us to look at genes beyond the ones that we normally look at.
    """
    df_cnv = df_calls[df_calls['variant_type'] == 'cnv'].copy()
    small_columns = ['runid', 'run_sample_id', 'call', 'variant_key_nt', 'percentage', 'gene']
    df_small_calls = df_cnv[small_columns]
    return df_small_calls


def override_calls(df_calls, df_reviewed, df_manifest):
    """
    Change manually adjudicated calls.

    Args;
        df_calls (pandas.DataFrame): all of the variant call information
        df_reviewed (pandas.DataFrame): the information for how to adjudicatethe calls
        df_manifest (pandas.DataFrame): the sample manifest. It should be the same one that was
            used to create ``df_calls``
    """
    df_new_calls = df_calls.copy()

    # Convert the wide format data frame to a tall format.
    # Each row is one call. One column has old call and one has the reviewed call
    mdf_reviewed = pd.melt(df_reviewed, id_vars=['patient_id', 'variant_key_nt'])
    mdf_reviewed['reviewed'] = mdf_reviewed['variable'].str.endswith('reviewed').\
        replace({True: 'reviewed', False: 'original'})
    mdf_reviewed['assay'] = mdf_reviewed['variable'].str.extract('(mdl|cdx|ldt)')
    df_tall = pd.pivot_table(mdf_reviewed, index=['patient_id', 'variant_key_nt', 'assay'],
                   columns='reviewed', values='value', aggfunc=lambda x: x).reset_index()
    df_tall = pd.merge(df_tall, df_manifest[['patient_id', 'assay', 'run_sample_id']], how='left')
    mdf_changed = df_tall[df_tall['original'] != df_tall['reviewed']]

    # Loop through each changed call
    for _, call_row in mdf_changed.iterrows():
        print(call_row['variant_key_nt'])
        assert call_row['reviewed'] in [0, 1, 2, 'remove']
        is_call = ((df_new_calls['run_sample_id'] == call_row['run_sample_id']) &
                   (df_new_calls['variant_key_nt'] == call_row['variant_key_nt']))
        n_calls = sum(is_call)
        assert (n_calls == 1) | ((call_row['original'] == 0) and (n_calls == 0))
        # Delete calls entirely if they were marked as removed or downgraded
        if call_row['reviewed'] in [0, 'remove']:
            df_new_calls = df_new_calls.loc[~is_call]
        # Update calls where there was a discrepancy in reported versus detected
        elif call_row['reviewed'] in [1, 2] and call_row['original'] in [1, 2]:
            df_new_calls.loc[is_call, 'pa_call'] = call_row['reviewed']
        # Create a call if it wasn't initially detected
        elif call_row['reviewed'] in [1, 2] and call_row['original'] == 0:
            template_call = df_calls[
                df_calls['variant_key_nt'] == call_row['variant_key_nt']]
            if template_call.empty:
                raise Exception('There is no template call for creating this call')
            template_call = template_call.iloc[0].copy()
            template_call['run_sample_id'] = call_row['run_sample_id']
            template_call['pa_call'] = call_row['reviewed']
            for none_column in ['call', 'percentage', 'variant_comment', 'somatic_call',
                                'fams', 'mut_fams']:
                template_call[none_column] = None
            sample_info = df_manifest[
                df_manifest['run_sample_id'] == call_row['run_sample_id']].iloc[0]
            for sample_column in ['unique_id', 'data_location', 'runid']:
                template_call[sample_column] = sample_info[sample_column]
            df_new_calls = df_new_calls.append(template_call)
        else:
            assert False
    return df_new_calls


def annotate_homopolymers(df_calls):
    """
    The protocol calls for some pretty specific variant classes. This function adds columns to
    annotate variants that have a 5-mer homopolymer within 5 bases of the edge of an indel.
    """
    df_indels = df_calls[df_calls['variant_type'] == 'indel'].copy()
    df_position = df_indels[['chrom', 'position', 'del_size', 'variant_type']].drop_duplicates()
    df_position['stop'] = df_position['position'].astype(int) + \
        df_position['del_size'].astype(int)

    for i_row, row in df_position.iterrows():
        chrom = row['chrom']

        # Record if the 5 bases before and after each edge of the the indel have a homopolymer.
        for edge in ['position', 'stop']:
            flank = hg19.fetch(f'chr{chrom}', int(row[edge]) - 4, int(row[edge]) + 6).upper()
            max_homopolymer = 1
            homopolymer_length = 1
            previous_letter = 'X'
            for letter in flank:
                if letter == previous_letter:
                    homopolymer_length += 1
                else:
                    homopolymer_length = 1
                previous_letter = letter
                if homopolymer_length > max_homopolymer:
                    max_homopolymer = homopolymer_length

        df_position.loc[i_row, 'homopolymer'] = max_homopolymer > 4

    df_position = df_position.drop('stop', axis=1)
    is_homopolymer = pd.merge(df_calls, df_position,
                              how='left')['homopolymer'].fillna(False).values
    return is_homopolymer


def annotate_long_indels(df_calls):
    """ Annotate which calls are a long indel """
    del_size = df_calls['del_size']
    ins_size = df_calls['ins_str'].str.len()
    return (del_size > 30) | (ins_size > 30)


def get_reportable_status(df_calls, key_column='variant_key_nt'):
    """
    Get a table indicating which variants are or are not in the reportable range.  This should be
    refactored to a function that is applied to only one bed file. It's too confusing to understand
    as written.

    .. code::

                                  GH2.10.1  GH2.10  GH2.11  MDACC_v1.0
        variant_key_nt
        ALK.2.29416572.T>C            True    True    True        True
        BRCA2.13.32915005.G>C         True    True    True        True
        BRCA2.13.32929387.T>C         True    True    True        True
        EML4.2.42524781.G>A          False   False   False       False
        FGFR2.10.123242026.T>C       False   False   False       False

    Args:
        df_calls (pandas.DataFrame): output of :py:func:`get_study_calls`

    Returns:
        pandas.DataFrame: The index is the `key_column` string and the columns are the panel
            names
    """
    panels = df_calls['panel'].unique()
    sftp_client = start_backend('sftp')

    # Separate the BED files by variant type
    vtype_beds = {panel: dict() for panel in BED_FILE_PATHS}
    for panel, file_path in BED_FILE_PATHS.items():
        with sftp_client.open(file_path, 'r') as file_pointer:
            bed = BedTool(file_pointer)
            vtype_intervals = {'snv': [], 'indel': [], 'noncoding': []}
            for interval in bed.features():
                info = interval.fields[5]
                info_dict = dict(item.split('=') for item in info.split(';'))
                reportable_vtypes = info_dict.get('REPORT', '').split(',')
                if 'snv' in reportable_vtypes:
                    vtype_intervals['snv'].append(interval)
                if 'indel' in reportable_vtypes:
                    vtype_intervals['indel'].append(interval)
                if 'noncoding' in info_dict.get('CAT', ''):
                    vtype_intervals['indel'].append(interval)
            for vtype, intervals in vtype_intervals.items():
                vtype_beds[panel][vtype] = BedTool(intervals)

    # Also get the blacklist bed files
    vtype_blacklist_beds = {panel: dict() for panel in BLACKLIST_BED_FILE_PATHS}
    for panel, file_path in BLACKLIST_BED_FILE_PATHS.items():
        with sftp_client.open(file_path, 'r') as file_pointer:
            bed = BedTool(file_pointer).remove_invalid()
            vtype_intervals = {'snv': [], 'indel': []}
            for interval in bed.features():
                snv_indel_note = interval.fields[3]
                if snv_indel_note == 'BLACKLIST-SNV':
                    vtype_intervals['snv'].append(interval)
                elif snv_indel_note == 'BLACKLIST-INDEL':
                    vtype_intervals['indel'].append(interval)
                else:
                    raise Exception('Unknown blacklist type')
            for vtype, intervals in vtype_intervals.items():
                vtype_blacklist_beds[panel][vtype] = BedTool(intervals)

    # For each variant type, create a BED file from all of the called variants. Then intersect it
    # with the reportable range BED file
    df_panels = df_calls[['chrom', 'position', 'del_size', key_column, 'variant_type',
                          'splice_effect', 'is_cs']].drop_duplicates().copy()
    for panel, panel_beds in vtype_beds.items():
        df_panels[panel] = False
        for vtype, vtype_bed in panel_beds.items():
            # Make a data frame and BED file of the relevant calls
            if vtype == 'snv':
                df_calls_vtype = df_panels[df_panels['variant_type'] == 'snv']
            elif vtype == 'indel':
                df_calls_vtype = df_panels[df_panels['variant_type'] == 'indel']
            elif vtype == 'noncoding':
                is_splice_site = df_panels['splice_effect'].str.startswith('splice').fillna(False)
                df_calls_vtype = df_panels[is_splice_site]
            # The positions from the call files are incremented by one to account of one-off
            # indexing difference between bed files and calls
            df_bed = pd.DataFrame({'chrom': df_calls_vtype['chrom'],
                                   'start': df_calls_vtype['position'].astype(int) + 1,
                                   'stop': df_calls_vtype['position'].astype(int) +
                                           df_calls_vtype['del_size'].astype(int) + 1,
                                   key_column: df_calls_vtype[key_column]})
            df_bed = df_bed[['chrom', 'start', 'stop', key_column]].\
                drop_duplicates().reset_index(drop=True)
            call_bed = BedTool.from_dataframe(df_bed)

            # Don't forget the blacklist
            # The "+" operator for BedTool objects returns an intersection of variants
            if vtype != 'noncoding':
                blacklist_bed = vtype_blacklist_beds[panel][vtype]
                reportable_calls = call_bed + vtype_bed - blacklist_bed
            else:
                reportable_calls = call_bed + vtype_bed

            try:
                reportable_variant_keys = reportable_calls.to_dataframe()['name'].tolist()
            except pd.errors.EmptyDataError:
                reportable_variant_keys = []
            df_panels.loc[df_panels[key_column].\
                              isin(reportable_variant_keys), panel] = True
            # Let's assume that fusions and CNVs are reportable if they are clinically significant
            is_fusion_cnv = df_panels['variant_type'].isin(['fusion', 'cnv'])
            df_panels.loc[is_fusion_cnv & df_panels['is_cs'], panel] = True
    df_panels = df_panels[[key_column] + list(panels)].drop_duplicates()
    df_reportable = df_panels.set_index(key_column)
    return df_reportable
