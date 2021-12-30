# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time

step_sizes = np.logspace(12, 14, 15)[::-1]
df_sim_data = pd.DataFrame(columns=('step_size', 'time', 'kappa', 'kappa_err'))
path = 'compact_sources/'
prop_module = 'BP'

for i, step_size in enumerate(step_sizes):
    crp = CRPropa(step_size = step_size, traj_max = 10**15, path = path, prop_module = prop_module)
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa, kappa_err = crp.analyze(step_size)

    df_sim_data.loc[i] = [step_size, time_needed, kappa, kappa_err]


# save simulation result
file_name_results = path + 'data/crp_sim_data_'+prop_module+'.pkl'
df_sim_data.to_pickle(file_name_results)
print(df_sim_data)
    