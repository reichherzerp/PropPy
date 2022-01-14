# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time
pc = 3.086*10**16

#step_sizes = np.logspace(5, 8, 7)[::-1]*pc
step_sizes = [10**4*pc]
df_sim_data = pd.DataFrame(columns=('step_size', 'time', 'kappa', 'kappa_err'))
path = 'uhecrs/'
prop_module = 'BP'
kappa_theory = 5.6*10**32 # [m^2/s]

# save simulation result
file_name_results = path + 'data/crp_sim_data_'+prop_module+'_test.pkl'

for i, step_size in enumerate(step_sizes[:10]):
    crp = CRPropa(energy = 1e19, brms = 10**(-9), l_max = 5*10**6*pc, l_min = 5*10**4*pc, step_size = step_size, traj_max = 10**11*pc, path = path, prop_module = prop_module, kappa = kappa_theory)
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa, kappa_err = crp.analyze(step_size)

    df_sim_data.loc[i] = [step_size, time_needed, kappa, kappa_err]
    df_sim_data.to_pickle(file_name_results) # save intermediate results
 

df_sim_data.to_pickle(file_name_results) # save final result
print(df_sim_data)
    