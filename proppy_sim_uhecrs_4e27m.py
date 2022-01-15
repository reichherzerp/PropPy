import numpy as np
import pandas as pd
import time
import rwpropa as pp
from pathlib import Path

path_figs = 'comparison/uhecrs_4e27m/figures'
path_data = 'comparison/uhecrs_4e27m/data'
path_data_raw = 'comparison/uhecrs_4e27m/data/raw_data'
Path(path_figs).mkdir(parents=True, exist_ok=True)
Path(path_data).mkdir(parents=True, exist_ok=True)
Path(path_data_raw).mkdir(parents=True, exist_ok=True)
file_name_results = path_data+'/proppy_sim_data.pkl'

pc = 3.086*10**16 # [m]

l_c = 10**6 # [pc]
energy = 10**19 # [eV]
r_g = 1.5*10**7 # [pc]
kappa_theory = 1.14 # [pc^2/s]
lambda_theory =  3.51*10**8 # [pc]
traj_max = 4e27/pc/10 # [pc]


step_sizes = ((np.logspace(20, 26, 19))/pc)


df_sim_data = pd.DataFrame(columns=('step_size', 'time', 'kappa', 'kappa_err'))

# start with large step sizes first, as these simulations are faster
step_sized_reversed = np.insert(step_sizes[::-1], 0, step_sizes[-1], axis=0)
for i, step_size in enumerate(step_sized_reversed):
    print('______________________________________________________________')
    print('step_size: ', step_size)
    nr_steps = int(traj_max/step_size)
    sim = pp.IsotropicSimulation(nr_particles = 10**2, energy = energy, step_size = step_size, nr_steps = nr_steps, diffusion_coefficient_para = kappa_theory, nr_obs_steps = 100, constants = pp.Constants(unit_distance = 2))
    file_name_raw = path_data_raw+'/proppy_stepsize_'+str(step_size/10**11)
    start_time = time.process_time()
    sim.simulate(file_name_raw)
    time_needed = time.process_time() - start_time
    print('time needed: ', time_needed, 's')
    df = pd.read_pickle(file_name_raw+'.pkl')
    sta = pp.Statistics(df)
    df_kappas = sta.get_diffusion_coefficients(speed_light=9.716*10**(-9)) # c in pc/s
    if i == 0:
        print('______________________________________________________________')
        print('Finished setup test - starting now with the convergence test!')
        continue
    df_sim_data.loc[i-1] = [step_size, time_needed, np.mean(df_kappas['kappa'][-10:]), np.std(df_kappas['kappa'][-10:])]
    df_sim_data.to_pickle(file_name_results)
    file_name = path_data+'/sim_result_proppy_stepsize_'
    np.save(file_name+str(step_size/10**11)+'_l.npy', np.array(df_kappas['l']))
    np.save(file_name+str(step_size/10**11)+'_kappa.npy', np.array(df_kappas['kappa']))
