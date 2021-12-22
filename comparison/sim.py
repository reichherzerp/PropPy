from crpropa_sim import *
import numpy as np

step_sizes = np.logspace(10**10, 10**12, 4)
crpropa_sim(10**11)
crpropa_sim(10**12)

analyze_agn(10**11, 'sim_result_ana_')
analyze_agn(10**12, 'sim_result_ana_')