# First, activate CRPropa environment and then run this current script
# Use the branch https://github.com/reichherzerp/CRPropa3/tree/LogObserverTimeEvolution
# to have access to the log ObserverTimeEvolution needed for the script!

from crpropa_sim import *
import numpy as np
import time
from pathlib import Path

path = 'compact_sources_4e27m_seed_analysis/'
path_figs = 'compact_sources_4e27m_seed_analysis/figures'
path_data = 'compact_sources_4e27m_seed_analysis/data'
path_data_raw = 'compact_sources_4e27m_seed_analysis/data/raw_data'
Path(path_figs).mkdir(parents=True, exist_ok=True)
Path(path_data).mkdir(parents=True, exist_ok=True)
Path(path_data_raw).mkdir(parents=True, exist_ok=True)

step_size = 3e21 # [m]
nr_seeds = 100

prop_module = 'BP'

# save simulation result
file_name_results = path + 'data/crp_sim_data_'+prop_module+'_seed_'+'.pkl'

df_sim_data = pd.DataFrame(columns=('seed', 'time', 'kappa', 'kappa_err'))

for i, seed in enumerate(range(nr_seeds)):
    crp = CRPropa(energy=1e19, brms=10**(-9), step_size = step_size, l_min = 5*3*10**20, l_max = 5*3e22, traj_max = 4e27, path = path, prop_module = prop_module, turbulence_method = 'PW', seed_study = True, random_seed = seed)
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa, kappa_err = crp.analyze(step_size)

    df_sim_data.loc[i] = [seed, time_needed, kappa, kappa_err]
    df_sim_data.to_pickle(file_name_results) # save intermediate results
 
print(df_sim_data)