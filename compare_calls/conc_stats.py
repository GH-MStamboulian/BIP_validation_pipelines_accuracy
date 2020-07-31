"""
Various functions to calculate metrics from contingency tables.
"""
import numpy as np
from numba import njit, prange
from statsmodels.stats.proportion import proportion_confint


def calc_agreement_stats(conc_vector, n_tests, n_boot, prior_freq):
    """
    Calculate several sorts of agreement metrics

    Args:
        conc_vector (numpy.ndarray): concordance vector
        n_tests (int): number of tests being compared
        n_boot (int): number of iterations in bootstrap calculations
    """
    overall = conc_vector_to_overall(conc_vector, n_tests)
    raw_ppa_npa_opa = calc_ppa_npa_opa(overall, n_tests)
    raw_bootstrap_values = run_bootstrap(conc_vector, n_tests, n_boot=n_boot)
    raw_bootstrap_cis = calc_bootstrap_cis(raw_bootstrap_values)
    out = {
        'overall': overall,
        'raw_ppa': raw_ppa_npa_opa[0],
        'raw_npa': raw_ppa_npa_opa[1],
        'raw_opa': raw_ppa_npa_opa[2],
        'raw_bootstrap_values': raw_bootstrap_values,
        'raw_ppa_llci_bootstrap': raw_bootstrap_cis[0],
        'raw_ppa_ulci_bootstrap': raw_bootstrap_cis[1],
        'raw_npa_llci_bootstrap': raw_bootstrap_cis[2],
        'raw_npa_ulci_bootstrap': raw_bootstrap_cis[3],
        'raw_opa_llci_bootstrap': raw_bootstrap_cis[4],
        'raw_opa_ulci_bootstrap': raw_bootstrap_cis[5],
        'n_boot': n_boot
    }
    out.update(calc_overall_stats(overall, n_tests))
    if n_tests == 3 and prior_freq is not None:
        cond_ppa_opa_npa = calc_ppa_npa_opa(overall, n_tests, prior_freq)
        cond_bootstrap_values = run_bootstrap(
            conc_vector, n_tests, n_boot=n_boot, prior_freq=prior_freq)
        cond_bootstrap_cis = calc_bootstrap_cis(cond_bootstrap_values)
        out.update({
            'cond_ppa': cond_ppa_opa_npa[0],
            'cond_npa': cond_ppa_opa_npa[1],
            'cond_opa': cond_ppa_opa_npa[2],
            'cond_bootstrap_values': cond_bootstrap_values,
            'cond_ppa_llci_bootstrap': cond_bootstrap_cis[0],
            'cond_ppa_ulci_bootstrap': cond_bootstrap_cis[1],
            'cond_npa_llci_bootstrap': cond_bootstrap_cis[2],
            'cond_npa_ulci_bootstrap': cond_bootstrap_cis[3],
            'cond_opa_llci_bootstrap': cond_bootstrap_cis[4],
            'cond_opa_ulci_bootstrap': cond_bootstrap_cis[5],
            'theta': calc_theta(overall, prior_freq)
        })

    # Now calculate the raw statistics
    if n_tests == 3:
        # Act as if the 3-test comparison is a 2-test comparison for the raw stats
        a0, b0, c0, d0, a1, b1, c1, d1 = overall
        a, b, c, d = a0 + a1, b0 + b1, c0 + c1, d0 + d1
        overall = [a, b, c, d]

    binomial_cis = calc_binomial_cis(overall)
    out.update({
        'raw_ppa_llci_binomial': binomial_cis[0],
        'raw_ppa_ulci_binomial': binomial_cis[1],
        'raw_npa_llci_binomial': binomial_cis[2],
        'raw_npa_ulci_binomial': binomial_cis[3],
        'raw_opa_llci_binomial': binomial_cis[4],
        'raw_opa_ulci_binomial': binomial_cis[5],
    })
    # What if we switched the role of reference and test?
    a, b, c, d = overall
    reverse_overall = [a, c, b, d]
    raw_ppv_npv_opv = calc_ppa_npa_opa(reverse_overall, n_tests=2)
    reverse_binomial_cis = calc_binomial_cis(reverse_overall)
    out.update({
        'raw_ppv': raw_ppv_npv_opv[0],
        'raw_npv': raw_ppv_npv_opv[1],
        'raw_opv': raw_ppv_npv_opv[2],
        'raw_ppv_llci_binomial': reverse_binomial_cis[0],
        'raw_ppv_ulci_binomial': reverse_binomial_cis[1],
        'raw_npv_llci_binomial': reverse_binomial_cis[2],
        'raw_npv_ulci_binomial': reverse_binomial_cis[3],
        'raw_opv_llci_binomial': reverse_binomial_cis[4],
        'raw_opv_ulci_binomial': reverse_binomial_cis[5],
    })
    return out


def calc_overall_stats(overall, n_tests):
    """
    Additional stats performed on the overall concordance counts.

    Args:
        overall (numpy.ndarray): from :py:func:`conc_vector_to_overall`
        n_tests (int): the number of tests being compared

    Returns:
        dict: the stats
    """
    if n_tests == 2:
        a, b, c, d = overall  # noqa
        ref_pos = a + b
        ref_neg = c + d
    elif n_tests == 3:
        a0, b0, c0, d0, a1, b1, c1, d1 = overall  # noqa
        ref_pos = a0 + b0 + a1 + b1
        ref_neg = c0 + d0 + c1 + d1
    out = {
        'ref_pos': ref_pos,
        'ref_neg': ref_neg
    }
    return out


def calc_binomial_cis(overall):
    """
    Calculate clopper-pearson confidence intervals given overall counts for 2-test comparsion
    """
    a, b, c, d = overall  # noqa
    ppa_llci, ppa_ulci = proportion_confint(a, a + b, method='beta')
    npa_llci, npa_ulci = proportion_confint(d, c + d, method='beta')
    opa_llci, opa_ulci = proportion_confint(a + d, a + b + c + d, method='beta')
    return [ppa_llci, ppa_ulci, npa_llci, npa_ulci, opa_llci, opa_ulci]


@njit
def conc_vector_to_overall(conc_vector, n_tests, random_patients=False):
    """
    Convert a flattened concordance vector that contains concordance counts for many patients into
    a smaller vector summed over all patients.

    Args:
        conc_vector (numpy.ndarray): concordance vector
        n_tests (int): number of tests being compared
        random_patients (bool): if ``True``, then randomly sample patients with replacement. This
            is useful as a base function for bootstrapping.

    Returns:
        numpy.ndarray: a vector of length 2**n_tests
    """
    n_cells = int(2**n_tests)  # Number of cells in contingency table
    n_patients = int(len(conc_vector) / n_cells)
    overall = [0] * n_cells

    if random_patients:
        patient_list = np.random.randint(0, n_patients - 1, n_patients)
    else:
        patient_list = np.arange(n_patients)
    for i_patient in patient_list:
        start_point = i_patient * n_cells
        for i_cell in range(n_cells):
            overall[i_cell] += conc_vector[start_point + i_cell]
    return overall


@njit
def safe_divide(numerator, denominator):
    return numerator / denominator if denominator > 0 else np.nan


def calc_theta(overall, prior_freq):
    """
    Just for calculating an intermediate statistic for display in a table. When it needs to be
    used, it's actually calculated in :py:func:`calc_ppa_npa_opa`.
    """
    a0, b0, c0, d0, a1, b1, c1, d1 = overall
    n0 = a0 + b0 + c0 + d0
    n1 = a1 + b1 + c1 + d1
    theta = safe_divide(n1 * (1-prior_freq), n0 * prior_freq)
    return theta


@njit(parallel=True)
def calc_ppa_npa_opa(overall, n_tests, prior_freq=None):
    """
    Calculate PPA, NPA, and OPA given the overall number of concordance types.

    Args:
        overall (numpy.ndarray): concordance types summed overall patients or output of
            :py:func:`conc_vector_to_overall`
        n_tests (int): number of tests being compared

    Returns:
        list: three floats. In order, the PPA, NPA, and OPA
    """
    if n_tests == 2:
        a, b, c, d = overall
        ppa = safe_divide(a, a + b)
        npa = safe_divide(d, c + d)
        opa = safe_divide(a + d, a + b + c + d)
    elif n_tests == 3:
        a0, b0, c0, d0, a1, b1, c1, d1 = overall
        n0 = a0 + b0 + c0 + d0
        n1 = a1 + b1 + c1 + d1
        if prior_freq is not None:
            theta = safe_divide(n1 * (1-prior_freq), n0 * prior_freq)
        else:
            theta = 1.0
        ppa = safe_divide(a0 * theta + a1, (a0 + b0) * theta + a1 + b1)
        npa = safe_divide(d0 * theta + d1, (c0 + d0) * theta + c1 + d1)
        opa = safe_divide((a0 + d0) * theta + a1 + d1, n0 * theta + n1)
    return [ppa, npa, opa]


@njit(parallel=True)
def run_bootstrap(conc_vector, n_tests, n_boot, prior_freq=None):
    """
    Run bootstrap simulations to get PPA, NPA, and OPA values.

    Args:
        conc_vector (numpy.ndarray): from :py:func:`conc_vector_to_overall`
        n_tests (int): number of tests being compared
        n_boot (int): number of bootstrap iterations
        prior_freq (float): prior frequency passed to :py:func:`calc_ppa_npa_opa`

    Returns:
        numpy.ndarray: the values from all of the bootstrap iterations. For example, for 2
        iterations, the vector would be [PPA1, NPA1, OPA1, PPA2, NPA2, OPA2].
    """
    n_metrics = 3
    values = np.zeros(n_boot * n_metrics)
    for i_boot in prange(n_boot):
        overall = conc_vector_to_overall(conc_vector, n_tests, random_patients=True)
        ppa_npa_opa = calc_ppa_npa_opa(overall, n_tests, prior_freq)
        for i_metric in range(n_metrics):
            values[i_boot * n_metrics + i_metric] = ppa_npa_opa[i_metric]
    return values


@njit
def calc_bootstrap_cis(bootstrap_values):
    """
    Calculate the confidence intervals via bootstrap sampling.

    Args:
        bootstrap_values (numpy.ndarray): from :py:func:`conc_vector_to_overall`

    Returns:
        list: six floats. In order, LLCI and ULCI for the PPA, then NPA, and then OPA
    """
    n_metrics = 3
    cis = [0.0] * n_metrics * 2
    for i_metric in range(n_metrics):
        metric_values = bootstrap_values[i_metric::3]
        metric_cis = np.nanpercentile(metric_values, [2.5, 97.5])
        cis[i_metric * 2] = metric_cis[0]  # LLCI
        cis[i_metric * 2 + 1] = metric_cis[1]  # ULCI
    return cis
