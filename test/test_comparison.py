import numpy as np
import proppy.comparison as comparison


def test_comparison():
    kappa_theory = 10**28
    lambda_theory = 10**23
    step_sizes = [10**21, 10**22]
    l_c = 10**19
    r_g = 10**19
    path_data_raw = ''
    path_data = ''
    path_figs = ''
    comp = comparison.Comparison(kappa_theory, lambda_theory, step_sizes, l_c, r_g, path_data_raw, path_data, path_figs, proppy_unit = 'm')

    comp.load_sim_data()