# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np


step_sizes = np.logspace(10, 13, 8)
for step_size in step_sizes:
    crp = CRPropa(step_size = step_size)
    crp.sim()
    crp.analyze(step_size, 'data/sim_result_crp_')