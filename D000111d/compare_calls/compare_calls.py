"""
Call comparisons
"""
import os
import numpy as np
import pandas as pd
from D000111d.settings import DATA_DIR, BED_FILE_PATHS
from pybedtools import BedTool


def get_call_comparisons(df_manifest, df_calls):
    """
    There are three sample collections.

    * Sample Collection 1: MDL Consecutive Samples
    * Sample Collection 2: MDL Enriched Positives
    * Sample Collection 3: GH Enrriched Positives

    This function gathers all of the manifest and call data and creates the comparison of calls for
    all the collections.

    patient_id          variant_key_nt  collection  is_cs ....
      A0133607      ALK.2.29416366.G>C           2  False
      A0133607      ALK.2.29416481.T>C           2  False

    clinical_class variant_type  cdx  mdl  ldt conc_type
                            snv    1    1   -1         d
                            snv    1    1   -1         d

    Args:
        df_manifest (pandas.DataFrame): the sample manifest
        df_calls (pandas.DataFrame): the annotated call information

    Returns:
        pandas.DataFrame: one row per patient-variant combination. The ``conc_type`` column
            indicates in which cell of a contingency table the comparison belongs. The call values
            in each assay column are the ``pa_call`` values (0=not detected, 1=detected,
            2=reported).
    """
    # Compare calls that are reportable in all assays
    df_12 = compare_reportable_calls(df_manifest[df_manifest['collection'].isin([1, 2])],
                                     df_calls[df_calls['collection'].isin([1, 2])])
    df_3 = compare_reportable_calls(df_manifest[df_manifest['collection'] == 3],
                                    df_calls[df_calls['collection'] == 3])
    df_comps = pd.concat([df_12, df_3], sort=False)
    df_comps['ldt'] = df_comps['ldt'].fillna(-1).astype(int)

    # Label which part of the contingency table should be filled in
    df_12 = classify_concordance(df_comps[df_comps['collection'].isin([1, 2])].copy(),
                                 ref_assay='mdl', test_assay='cdx')
    df_3 = classify_concordance(df_comps[df_comps['collection'] == 3].copy(),
                                ref_assay='mdl', test_assay='cdx', select_assay='ldt')
    df_comps = pd.concat([df_12, df_3], sort=False)
    return df_comps


def get_call_comparisons_2test(df_manifest, df_calls, ref_assay='mdl', test_assay='cdx'):
    """ Copy :py:func:`get_call_comparisons` except that don't condition on LDT results for the
    third sample collection """
    # Figure out which variants are in the reportable range for each panel
    df_manifest_2 = df_manifest[df_manifest['assay'].isin([ref_assay, test_assay])]
    df_calls_2 = df_calls[df_calls['assay'].isin([ref_assay, test_assay])]

    # Compare calls that are reportable in all assays
    df_comps = compare_reportable_calls(df_manifest_2, df_calls_2)
    df_comps = classify_concordance(df_comps, ref_assay=ref_assay, test_assay=test_assay)
    return df_comps


def compare_reportable_calls(df_manifest,
                             df_calls,
                             key_column='variant_key_nt'):
    """
    Take call data that is in tall format and pivot it to have a boolean column for each assay
    indicating if the variant is called positive or not for that assay. The assumption is that
    ``df_calls`` will have an 'assay' column to make this possible.

    A pivoted row will only be included if it is reportable in all assays. The information about
    reportable variants is provided with ``df_reportable``.  The ``pa_call`` encodes if the call
    was detected (1), detected as somatic (2), or not detected (0).

    Args:
        df_manifest (pandas.DataFrame): must have three columns: 'patient_id', 'assay', and 'panel'
        df_calls (pandas.DataFrame): the output of
            :py:func:`av360.extract.extract_calls.get_catted_call_tables` with additional columns:
            'patient_id', 'panel', and 'assay'
        key_column (str): the name of the column to uniquely identify a variant

    Returns:
        pandas.DataFrame: each call on a column. All the columns are...

        * ``patient_id``
        * ``{key_column}``
        * ``is_cs`` (is a clinically significant variant)
        * ``variant_type``
        * then there is an integer column for each assay (0=not detected, 1=detected, 2=somatic
          call)
        * ``conc_type``: a string indicating where in a contingency table this comparison belongs
          (e.g., a, b, c, d, a0, b0, ..., c0, d0)
    """
    # Only look at calls for samples in the manifest
    df_calls = df_calls[df_calls['run_sample_id'].isin(df_manifest['run_sample_id'])]

    # Make sure that every sample has a call (to make sure the call file was actually parsed).
    # This works because we assume that every sample has at least one call (not necessarily the
    # case)
    # assert df_calls['run_sample_id'].nunique() == len(df_manifest)

    # Figure out which variants are in the reportable range for each panel
    panels = df_manifest['panel'].unique()
    panel_cols = ['panel_' + panel for panel in panels] + [key_column]
    df_reportable = df_calls[panel_cols].drop_duplicates().set_index(key_column)
    df_reportable.columns = [col.replace('panel_', '') for col in df_reportable.columns]

    # Only use the columns that are necessary
    index_columns = ['patient_id', key_column, 'is_cs', 'clinical_class', 'collection',
                     'homopolymer', 'long_indel', 'variant_type']
    call_columns = index_columns + ['pa_call', 'assay']
    df_calls = df_calls[call_columns].drop_duplicates()

    # Compare calls across the ``'assay'``
    df_comp = pd.pivot_table(df_calls, index=index_columns, columns='assay', values='pa_call').\
        fillna(0).astype(int).reset_index()

    # Get the panel for every combination of assay, patient, and variant key
    # patient_id                      variant_key_nt       ldt     cdx         mdl
    #   A0132472               BRCA1.17.41244004.G>A  GH2.10.1  GH2.11  MDACC_v1.0
    assays = df_calls['assay'].unique().tolist()
    df_panels = pd.pivot(df_manifest, index='patient_id', columns='assay', values='panel')
    df_left = df_comp[['patient_id', key_column]]
    df_right = df_panels.reset_index()[['patient_id'] + assays]
    df_call_panels = pd.merge(df_left, df_right, how='left')

    # Label each each call as reportable (or not reportable) for the panel that was run
    # patient_id                      variant_key_nt    ldt    cdx    mdl
    #   A0132472               BRCA1.17.41244004.G>A   True   True   True
    df_call_reportable = df_call_panels.copy()
    mdf_reportable = pd.melt(df_reportable.reset_index(), id_vars=key_column,
                             var_name='panel', value_name='reportable')
    for assay in assays:
        df_left = df_call_panels[[key_column, assay]].rename(columns={assay: 'panel'})
        df_call_reportable[assay] = pd.merge(df_left, mdf_reportable, how='left')['reportable']

    # Make sure that everything is reportbale or not and there are no missing values
    assert not df_call_reportable.isnull().any().any()

    # Only look at calls that are reportable in all ``assays``
    assert all(df_call_reportable['patient_id'] == df_comp['patient_id'])
    assert all(df_call_reportable[key_column] == df_comp[key_column])
    all_reportable = df_call_reportable[assays].all(axis=1)
    df_reportable_comp = df_comp[all_reportable].copy()
    return df_reportable_comp


def classify_concordance(df_comp, ref_assay, test_assay, select_assay=None):
    """
    When calls are compared, they go into different cells of a contingency table. This function
    takes the calls as input and outputs a string indicating which cell of the contingency table it
    should go into.

    This function assumes that 0=not detected, 1=detected, and 2=reported.

    Args:
        df_comp (pandas.DataFrame): one row for each variant call. The call for each assay are the
            columns, from :py:func:`transform.compare_calls.get_call_comparisons`
        ref_assay (str): name of reference assay
        test_assay (str): name of test assay
        select_assay (str): name of assay used to select samples. If ``None`` (the default), then
            it is assume there was no "select assay"

    Returns:
        pandas.DataFrame: ``df_comp`` with one additional column, "conc_type"
    """
    ref, test = df_comp[ref_assay], df_comp[test_assay]
    if select_assay is not None:
        select = df_comp[select_assay]
        df_conc = pd.DataFrame({
            'a0': (((ref == 1) & (test == 2) & (select == 0)) |
                   ((ref == 1) & (test == 2) & (select == 1)) |
                   ((ref == 2) & (test == 1) & (select == 0)) |
                   ((ref == 2) & (test == 1) & (select == 1)) |
                   ((ref == 2) & (test == 2) & (select == 0)) |
                   ((ref == 2) & (test == 2) & (select == 1))),
            'b0': (((ref == 2) & (test == 0) & (select == 0)) |
                   ((ref == 2) & (test == 0) & (select == 1))),
            'c0': (((ref == 0) & (test == 2) & (select == 0)) |
                   ((ref == 0) & (test == 2) & (select == 1))),
            'd0': (((ref == 0) & (test == 0) & (select == 0)) |
                   ((ref == 0) & (test == 0) & (select == 1)) |
                   ((ref == 0) & (test == 1) & (select == 0)) |
                   ((ref == 0) & (test == 1) & (select == 1)) |
                   ((ref == 1) & (test == 0) & (select == 0)) |
                   ((ref == 1) & (test == 0) & (select == 1)) |
                   ((ref == 1) & (test == 1) & (select == 0)) |
                   ((ref == 1) & (test == 1) & (select == 1))),
            'a1': (((ref == 1) & (test == 2) & (select == 2)) |
                   ((ref == 2) & (test == 1) & (select == 2)) |
                   ((ref == 2) & (test == 2) & (select == 2))),
            'b1': (((ref == 2) & (test == 0) & (select == 2))),
            'c1': (((ref == 0) & (test == 2) & (select == 2))),
            'd1': (((ref == 0) & (test == 0) & (select == 2)) |
                   ((ref == 0) & (test == 1) & (select == 2)) |
                   ((ref == 1) & (test == 0) & (select == 2)) |
                   ((ref == 1) & (test == 1) & (select == 2)))
        })
    else:
        df_conc = pd.DataFrame({
            'a': (((ref == 2) & (test == 1)) |
                  ((ref == 2) & (test == 2)) |
                  ((ref == 1) & (test == 2))),
            'b': (ref == 2) & (test == 0),
            'c': (ref == 0) & (test == 2),
            'd': (((ref == 0) & (test == 0)) |
                  ((ref == 0) & (test == 1)) |
                  ((ref == 1) & (test == 0)) |
                  ((ref == 1) & (test == 1)))
        })

    cell = pd.Series(None, index=df_comp.index)
    for col in df_conc:
        cell[df_conc[col]] = col
    df_comp['conc_type'] = cell
    return df_comp


def make_call_count_vector(df_comps, patient_ids, n_vars, n_tests):
    """
    Create a flattened that have the number of concordances and discordances.

    Args:
        df_comps (pandas.DataFrame): must have the following columns: "conc_type" and "patient_id"
        patient_ids (array-like): a vector of patient IDs. This also sets the order of the returned
            vector
        n_vars (int): the number of variants assumes to exist
        n_tests (int): the number of tests being compared

    Returns:
        numpy.ndarray: a vector of concordance counts
    """
    # Set up variables for making vectors
    if n_tests == 2:
        measure_vars = ['a', 'b', 'c', 'd']
        all_neg_cell = 'd'
    elif n_tests == 3:
        measure_vars = ['a0', 'b0', 'c0', 'd0', 'a1', 'b1', 'c1', 'd1']
        all_neg_cell = 'd0'
    else:
        raise Exception('Only 2 or 3 tests allowed')

    # Count up everything that isn't negative in all tests
    mdf_counts = df_comps.groupby(['patient_id', 'conc_type']).size().reset_index()
    df_counts = pd.pivot(mdf_counts, index='patient_id', columns='conc_type', values=0).\
        fillna(0).astype(int)

    # Make sure that there is a column for every concordance type
    for measure_var in measure_vars:
        if measure_var not in df_counts:
            df_counts[measure_var] = 0

    # If a patient has only calls that are negative in all assays, make sure it's not left out
    missing_ids = set(patient_ids).difference(df_counts.index)
    df_missing = pd.DataFrame(0, index=missing_ids, columns=measure_vars)
    df_counts = pd.concat([df_counts, df_missing])

    # Make sure the data frame is ordered the same as the given patient_id vector
    df_counts = df_counts.loc[patient_ids]

    # The number of concordant negatives is the total number of variants minus other stuff
    df_counts[all_neg_cell] = \
        n_vars - df_counts.drop(all_neg_cell, axis=1, errors='ignore').sum(axis=1)
    count_vector = np.ndarray.flatten(df_counts[measure_vars].values)
    return count_vector
