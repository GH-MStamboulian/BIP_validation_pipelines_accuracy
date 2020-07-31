"""
We want to compare different slices and groups of patients and variants. These functions do all of
that.
"""
from functools import reduce
import numpy as np
from joblib import Memory
from D000111d.settings import CACHE_DIR, PRIORS
from D000111d.compare_calls.compare_calls import make_call_count_vector
from D000111d.compare_calls.conc_stats import calc_agreement_stats, calc_bootstrap_cis
from D000111d.compare_calls.conc_stats import calc_ppa_npa_opa, calc_binomial_cis


N_BOOT = 100000
memory = Memory(CACHE_DIR, verbose=0)
c_calc_agreement_stats = memory.cache(calc_agreement_stats)


def compare_collection_groups(df_manifest, df_comps, priors=PRIORS):
    """
    Given the manifest and list of concordances for all variants, calculate various statistics. See
    :py:func:`comparee_calls.conc_stats.calc_agreement_stats` for more info about which sorts of
    stats are calculated.

    Args:
        df_manifest (pandas.DataFrame): the sample manifest
        df_comps (pandas.DataFrame): list of variant comparisons

    Returns:
        list: each slice of patients/variants are packaged into a dict of statistics
    """
   # Make call comparison vectors for all different slices of patients and variants
    groups = []

    for collection in [1, 2, 3]:
        group_manifest = df_manifest[df_manifest['collection'] == collection]
        if group_manifest.empty:
            continue
        coll_comps = df_comps[df_comps['collection'] == collection]
        if set(['a', 'b', 'c', 'd']).intersection(coll_comps['conc_type']):
            n_tests = 2
        elif set(['a0', 'b0', 'c0', 'd0', 'a1', 'b1', 'c1', 'd1']).\
            intersection(coll_comps['conc_type']):
            n_tests = 3

        # First, the variant categories
        for (is_cs, variant_type), n_vars in priors['category']['n_var'].items():
            cs_name = 'cs' if is_cs else 'panel'
            group_name = f'coll{collection}_{cs_name}_{variant_type}'
            group_comps = coll_comps[coll_comps['variant_type'] == variant_type]
            if is_cs:
                group_comps = group_comps[group_comps['is_cs']]
            conc_vector = make_call_count_vector(
                group_comps, group_manifest['patient_id'].unique(), n_vars, n_tests)
            prior_freq = priors['category']['prior'][(is_cs, variant_type)]
            var_per_sample = priors['category']['var_per_sample'][(is_cs, variant_type)]
            groups.append({
                'name': group_name,
                'class': f'{cs_name}_{variant_type}',
                'vector': conc_vector,
                'collection': collection,
                'n_tests': n_tests,
                'stats': c_calc_agreement_stats(conc_vector, n_tests, N_BOOT, prior_freq),
                'section': 'variant_category',
                'prior_freq': prior_freq,
                'n_vars': n_vars,
                'var_per_sample': var_per_sample
            })

        # Then, the clinical classes
        for clinical_class, n_vars in priors['clinical_class']['n_var'].items():
            group_name = f'coll{collection}_{clinical_class}'
            group_comps = coll_comps[coll_comps['clinical_class'] == clinical_class]
            conc_vector = make_call_count_vector(
                group_comps, group_manifest['patient_id'].unique(), n_vars, n_tests)
            prior_freq = priors['clinical_class']['prior'][clinical_class]
            var_per_sample = priors['clinical_class']['var_per_sample'][clinical_class]
            groups.append({
                'name': group_name,
                'class': clinical_class,
                'vector': conc_vector,
                'collection': collection,
                'n_tests': n_tests,
                'stats': c_calc_agreement_stats(conc_vector, n_tests, N_BOOT, prior_freq),
                'section': 'variant_class',
                'prior_freq': prior_freq,
                'n_vars': n_vars,
                'var_per_sample': var_per_sample
            })

        # Last, the special variant classes (long indels and homopolymers)
        for special_class, n_vars in priors['special']['n_var'].items():
            group_name = f'coll{collection}_{special_class}'
            # The call data was configured to have a boolean column named according to the
            # "special_class"
            group_comps = coll_comps[coll_comps[special_class]]
            conc_vector = make_call_count_vector(
                group_comps, group_manifest['patient_id'].unique(), n_vars, n_tests)
            prior_freq = priors['special']['prior'][special_class]
            var_per_sample = priors['special']['var_per_sample'][special_class]
            groups.append({
                'name': group_name,
                'class': special_class,
                'vector': conc_vector,
                'collection': collection,
                'n_tests': n_tests,
                'stats': c_calc_agreement_stats(conc_vector, n_tests, N_BOOT, prior_freq),
                'section': 'special',
                'prior_freq': prior_freq,
                'n_vars': n_vars,
                'var_per_sample': var_per_sample
            })

        # Last minute addition, also include the CNVs and fusions. One variant at a time.
        single_gene_groups = [
            ('cnv', 'MET'),
            ('cnv', 'ERBB2'),
            ('fusion', 'ALK'),
            ('fusion', 'RET'),
            ('fusion', 'ROS1'),
            ('fusion', 'NTRK1'),
        ]
        for variant_type, gene in single_gene_groups:
            n_vars = 1
            special_class = f'{variant_type}:{gene}'
            group_name = f'coll{collection}_{special_class}'
            group_comps = coll_comps[
                (coll_comps['variant_type'] == variant_type) &
                coll_comps['variant_key_nt'].str.contains(gene)
            ]
            conc_vector = make_call_count_vector(
                group_comps, group_manifest['patient_id'].unique(), n_vars, n_tests)
            # Let's just set dummy `prior_freq` and `var_per_sample` values assuming that we will
            # never view the conditioned results for the single-gene groups
            prior_freq = 0.01
            var_per_sample = 0.01
            groups.append({
                'name': group_name,
                'class': special_class,
                'vector': conc_vector,
                'collection': collection,
                'n_tests': n_tests,
                'stats': c_calc_agreement_stats(conc_vector, n_tests, N_BOOT, prior_freq),
                'section': 'single_gene',
                'prior_freq': prior_freq,
                'n_vars': n_vars,
                'var_per_sample': var_per_sample
            })

    return groups


def compare_combined_collections(groups):
    """
    Combine sample collections to get overall statistics for each variant class.

    Args:
        groups (list): output of :py:func:`compare_collection_groups`

    Returns:
        list: similar format as the input
    """
    combined_groups = []

    # For each class, calculate a combined PPA and NPA and make a new "comparison group"
    classes = {group['class'] for group in groups}
    for variant_class in classes:
        new_group_stats = dict()
        collection_groups = [group for group in groups if group['class'] == variant_class]
        group_section = list({group['section'] for group in collection_groups})[0]

        # Simply a 3-test comparison to a 2-test comparison
        overall = np.zeros(4, dtype=int)
        for group in collection_groups:
            group_stats = group['stats']
            n_tests = group['n_tests']
            if n_tests == 2:
                overall = overall + np.array(group_stats['overall'])
            elif n_tests == 3:
                overall = overall + np.array(group_stats['overall'])[0:4]
                overall = overall + np.array(group_stats['overall'])[4:8]
        new_group_stats = {'overall': overall}

        # For binomial PPA and NPA, we have to consider everything as only two tests
        # i.e., we have to use the simplification above
        raw_ppa_npa_opa = calc_ppa_npa_opa(overall, n_tests=2)
        binomial_cis = calc_binomial_cis(overall)
        new_group_stats.update({
            f'raw_ppa': raw_ppa_npa_opa[0],
            f'raw_npa': raw_ppa_npa_opa[1],
            f'raw_opa': raw_ppa_npa_opa[2],
            f'raw_ppa_llci_binomial': binomial_cis[0],
            f'raw_ppa_ulci_binomial': binomial_cis[1],
            f'raw_npa_llci_binomial': binomial_cis[2],
            f'raw_npa_ulci_binomial': binomial_cis[3],
            f'raw_opa_llci_binomial': binomial_cis[4],
            f'raw_opa_ulci_binomial': binomial_cis[5]
        })

        # We can add together the bootstrapping results from each group. This will also make sure
        # that the combined bootstrapping results are consistent with the bootstrapping results
        # from each collection
        for correction_type in ['raw', 'cond']:
            numerators = []
            denoms = []
            point_numerators = []
            point_denoms = []
            # Loop through each sample collection for the variant class
            for group in collection_groups:
                group_stats = group['stats']
                n_tests = group['n_tests']

                # For 2-test comparisons, we are not doing any conditioning because there can't be
                if n_tests == 2 and correction_type == 'cond':
                    bootstrap_values = group_stats[f'raw_bootstrap_values']
                    point_estimates = np.array([
                        group_stats['raw_ppa'], group_stats['raw_npa'], group_stats['raw_opa']])
                else:
                    bootstrap_values = group_stats[f'{correction_type}_bootstrap_values']
                    point_estimates = np.array([
                        group_stats[f'{correction_type}_ppa'],
                        group_stats[f'{correction_type}_npa'],
                        group_stats[f'{correction_type}_opa']])

                # The weights for PPA and NPA values are the number of reference calls
                ref_pos, ref_neg, n_boot = \
                    group_stats['ref_pos'], group_stats['ref_neg'], group_stats['n_boot']
                n_calls = np.array([ref_pos, ref_neg, ref_pos + ref_neg])

                # Assuming metrics are stored in order as PPA, NPA, and OPA
                weights = np.tile(n_calls, n_boot)
                weighted_bootstrap_values = weights * bootstrap_values
                # If there are any NaN values because of no calls, then set them to zero
                weighted_bootstrap_values[np.isnan(bootstrap_values) & (weights == 0)] = 0.0
                numerators.append(weighted_bootstrap_values)
                denoms.append(weights)

                # Don't forget to calculate the point estimates
                point_sum = point_estimates * n_calls
                point_sum[np.isnan(point_sum) & (n_calls == 0)] = 0.0
                point_numerators.append(point_sum)
                point_denoms.append(n_calls)

            # Sum across all collections
            if len(numerators) == len(collection_groups) and len(denoms) == len(collection_groups):
                # Bootstrap values
                numerator = reduce(lambda a, b: a + b, numerators)
                denom = reduce(lambda a, b: a + b, denoms)
                combined_bootstrap_values = numerator / denom
                combined_cis = calc_bootstrap_cis(combined_bootstrap_values)

                new_group_stats.update({
                    f'{correction_type}_ppa_llci_bootstrap': combined_cis[0],
                    f'{correction_type}_ppa_ulci_bootstrap': combined_cis[1],
                    f'{correction_type}_npa_llci_bootstrap': combined_cis[2],
                    f'{correction_type}_npa_ulci_bootstrap': combined_cis[3],
                    f'{correction_type}_opa_llci_bootstrap': combined_cis[4],
                    f'{correction_type}_opa_ulci_bootstrap': combined_cis[5]
                })
                # Point estimates only for conditional PPA and NPA
                if correction_type == 'cond':
                    point_numerator = reduce(lambda a, b: a + b, point_numerators)
                    point_denom = reduce(lambda a, b: a + b, point_denoms)
                    combined_ppa_npa_opa = point_numerator / point_denom
                    new_group_stats.update({
                        f'cond_ppa': combined_ppa_npa_opa[0],
                        f'cond_npa': combined_ppa_npa_opa[1],
                        f'cond_opa': combined_ppa_npa_opa[2],
                    })
        new_group = {
            'n_boot': n_boot,
            'collection': 4,
            'class': variant_class,
            'name': f'coll4_{variant_class}',
            'section': group_section,
            'stats': new_group_stats,
        }
        combined_groups.append(new_group)
    return combined_groups
