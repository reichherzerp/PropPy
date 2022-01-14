# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time
from pathlib import Path

path = 'uhecrs_3e28m_seed_analysis/'
path_figs = 'uhecrs_3e28m_seed_analysis/figures'
path_data = 'uhecrs_3e28m_seed_analysis/data'
path_data_raw = 'uhecrs_3e28m_seed_analysis/data/raw_data'
Path(path_figs).mkdir(parents=True, exist_ok=True)
Path(path_data).mkdir(parents=True, exist_ok=True)
Path(path_data_raw).mkdir(parents=True, exist_ok=True)

step_size = 3e21 # [m]
kappa_theory = 1.4*1e31 # [m^2/s]
nr_seeds = 8

prop_module = 'BP'

# save simulation result
file_name_results = path + 'data/crp_sim_data_'+prop_module+'_seed_analysis_30EeV.pkl'

df_sim_data = pd.DataFrame(columns=('seed', 'wm', 'time', 'kappa', 'kappa_err'))

nr_wave_modes = [10,10,100,100,1000,1000,10000,10000]
trajectory_lengths = [1e30, 1e30, 1e29, 1e29, 1e29, 1e29, 1e28, 1e28]
for i, seed in enumerate(range(nr_seeds)):
    wm = nr_wave_modes[i]
    crp = CRPropa(energy = 30*10**19, brms=10**(-9), step_size = step_size, l_min = 5*3e20, l_max = 5*3e22, traj_max = trajectory_lengths[i], path = path, nr_grid_points = 1024, prop_module = prop_module, kappa = kappa_theory, turbulence_method = 'PW', seed_study = True, random_seed = seed, n_wavemodes = wm)
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa, kappa_err = crp.analyze(step_size)

    df_sim_data.loc[i] = [seed, wm, time_needed, kappa, kappa_err]
    df_sim_data.to_pickle(file_name_results) # save intermediate results
 
print(df_sim_data)