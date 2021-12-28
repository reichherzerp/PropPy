# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time

step_sizes = np.logspace(11, 14, 15)[::-1]
df_sim_data = pd.DataFrame(columns=('step_size', 'time', 'kappa'))

for i, step_size in enumerate(step_sizes):
    crp = CRPropa(step_size = step_size, traj_max = 10**16)
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa = crp.analyze(step_size, 'data/sim_result_crp_')

    df_sim_data.loc[i] = [step_size, time_needed, kappa]

# save simulation result
file_name_results = 'data/crp_sim_data.pkl'
df_sim_data.to_pickle(file_name_results)
print(df_sim_data)
    
