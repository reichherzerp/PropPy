# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time
from pathlib import Path

path = 'compact_sources_1e17m/'
path_figs = 'compact_sources_1e17m/figures'
path_data = 'compact_sources_1e17m/data'
path_data_raw = 'compact_sources_1e17m/data/raw_data'
Path(path_figs).mkdir(parents=True, exist_ok=True)
Path(path_data).mkdir(parents=True, exist_ok=True)
Path(path_data_raw).mkdir(parents=True, exist_ok=True)

step_sizes = np.logspace(9, 15, 19)
kappa_theory = 1.59*10**23 # [m^2/s]

prop_module = 'BP'
# save simulation result
file_name_results = path + 'data/crp_sim_data_'+prop_module+'.pkl'

df_sim_data = pd.DataFrame(columns=('step_size', 'time', 'kappa', 'kappa_err'))

for i, step_size in enumerate(step_sizes[::-1]):
    crp = CRPropa(step_size = step_size, l_min = 5*10**9, nr_grid_points = 1024, traj_max = 10**17, path = path, prop_module = prop_module, kappa = kappa_theory, turbulence_method = 'grid')
    start_time = time.process_time()
    crp.sim()
    time_needed = time.process_time() - start_time
    
    kappa, kappa_err = crp.analyze(step_size)

    df_sim_data.loc[i] = [step_size, time_needed, kappa, kappa_err]
    df_sim_data.to_pickle(file_name_results) # save intermediate results
 
print(df_sim_data)