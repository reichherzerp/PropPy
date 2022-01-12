# First, activate CRPropa environment and then run this current script

from crpropa_sim import *
import numpy as np
import time
from pathlib import Path


path_figs = 'compact_sources_1e17m/figures'
path_data = 'compact_sources_1e17m/data'
path_data_raw = 'compact_sources_1e17m/data/raw_data'
Path(path_figs).mkdir(parents=True, exist_ok=True)
Path(path_data).mkdir(parents=True, exist_ok=True)
Path(path_data_raw).mkdir(parents=True, exist_ok=True)
# save simulation result
path = 'compact_sources_1e17m/'

step_sizes = [(np.logspace(9, 15, 19))[8]]
df_sim_data = pd.DataFrame(columns=('step_size', 'time', 'kappa', 'kappa_err'))

kappa_theory = 3.34*10**22 # [m^2/s]
traj_max = 10**17 # [m]

simulation_setups = [
    {
        'prop_module': 'SDE',
        'turbulence_method': '',
        'nr_grid_points': 0,
        'nr_modes': 0
    },
    {
        'prop_module': 'BP',
        'turbulence_method': 'PW',
        'nr_grid_points': 0,
        'nr_modes': 1000
    },
    {
        'prop_module': 'CK',
        'turbulence_method': 'PW',
        'nr_grid_points': 0,
        'nr_modes': 1000
    },
    {
        'prop_module': 'BP',
        'turbulence_method': 'grid',
        'nr_grid_points': 256,
        'nr_modes': 0,
    },
]

def simulate(simulation_setup):
    file_name_results = path + 'data/crp_sim_data_'+simulation_setup['prop_module']+'_'+simulation_setup['turbulence_method']+'.pkl'
    for i, step_size in enumerate(step_sizes):
        crp = CRPropa(step_size = step_size, traj_max = traj_max, nr_grid_points = simulation_setup['nr_grid_points'], n_wavemodes = simulation_setup['nr_modes'], l_min = 5*10**9, path = path, prop_module = simulation_setup['prop_module'], kappa = kappa_theory, turbulence_method = simulation_setup['turbulence_method'])
        start_time = time.process_time()
        crp.sim()
        time_needed = time.process_time() - start_time
        
        kappa, kappa_err = crp.analyze(step_size)

        df_sim_data.loc[i] = [step_size, time_needed, kappa, kappa_err]
        df_sim_data.to_pickle(file_name_results) # save intermediate results
    
    print(df_sim_data)
 
for simulation_setup in [simulation_setups[0]]:
    simulate(simulation_setup)