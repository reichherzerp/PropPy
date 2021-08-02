import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from files.Simulation import Simulation
from files.Source import Source
from files.Observer import Observer
from files.Propagator import Propagator
from files.Observer import TimeEvolutionObserverLog
from plot.Trajectory import Trajectory
from plot.Statistics import Statistics

sim = Simulation()

nr_particles = 10**1
source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
delta_rho_div_phi = 1 # (delta_r_rho / delta_r_phi)
gyro_radius = 10**11 # [m]; 1pc

source = Source(gyro_radius, source_pos, nr_particles)
sim.add_source(source)

nr_steps = 10**4
step_size = 0.5*10**10 # [m]
mfp = np.array([3.75*10**13/4.0, 3.75*10**13/4.0, 7.2*10**13], dtype=np.float32)  # [m]

propagator = Propagator(nr_steps, step_size, mfp)
sim.add_propagator(propagator)

substeps = [False, False, True] # observe only steps (no substeps)
min_step = 1
max_step = nr_steps
nr_obs_steps = 100

observer = TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)
sim.add_observer(observer)

print('start sim')
sim.run_simulation()
print('finished sim')
sim.save_data('data')