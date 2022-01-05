# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time

step_sizes = np.logspace(10, 15, 15)[:5]
df_sim_data = pd.DataFrame(columns=('step_size', 'time', 'kappa', 'kappa_err'))
path = 'compact_sources_1e17m/'
prop_module = 'SDE'
kappa_theory = 1.59*10**23 # [m^2/s]

# save simulation result
file_name_results = path + 'data/crp_sim_data_'+prop_module+'.pkl'

for i, step_size in enumerate(step_sizes):
    crp = CRPropa(step_size = step_size, traj_max = 10**17, path = path, prop_module = prop_module, kappa = kappa_theory)
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa, kappa_err = crp.analyze(step_size)

    df_sim_data.loc[i] = [step_size, time_needed, kappa, kappa_err]
    df_sim_data.to_pickle(file_name_results) # save intermediate results
 

df_sim_data.to_pickle(file_name_results) # save final result
print(df_sim_data)
    