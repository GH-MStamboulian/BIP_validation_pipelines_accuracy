# first line: 9
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
