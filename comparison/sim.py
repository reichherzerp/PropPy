from crpropa_sim import *
import numpy as np

step_sizes = np.logspace(10, 13, 8)
for step_size in step_sizes:
    crpropa_sim(step_size = step_size)

for step_size in step_sizes:
    analyze_agn(step_size, 'data/sim_result_ana_')