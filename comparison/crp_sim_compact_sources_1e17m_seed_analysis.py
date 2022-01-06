# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time
from pathlib import Path

path = 'compact_sources_1e17m_seed_analysis/'
path_figs = 'compact_sources_1e17m_seed_analysis/figures'
path_data = 'compact_sources_1e17m_seed_analysis/data'
path_data_raw = 'compact_sources_1e17m_seed_analysis/data/raw_data'
Path(path_figs).mkdir(parents=True, exist_ok=True)
Path(path_data).mkdir(parents=True, exist_ok=True)
Path(path_data_raw).mkdir(parents=True, exist_ok=True)

step_size = 10**11 # [m]
kappa_theory = 1.59*10**23 # [m^2/s]
nr_seeds = 5

prop_module = 'BP'
# save simulation result
file_name_results = path + 'data/crp_sim_data_'+prop_module+'_seed_'+'.pkl'

df_sim_data = pd.DataFrame(columns=('seed', 'time', 'kappa', 'kappa_err'))

for i, seed in enumerate(range(nr_seeds)):
    crp = CRPropa(step_size = step_size, l_min = 5*10**9, traj_max = 10**16, path = path, prop_module = prop_module, kappa = kappa_theory, turbulence_method = 'PW', seed_study=True, seed= seed)
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa, kappa_err = crp.analyze(step_size)

    df_sim_data.loc[i] = [seed, time_needed, kappa, kappa_err]
    df_sim_data.to_pickle(file_name_results) # save intermediate results
 
print(df_sim_data)