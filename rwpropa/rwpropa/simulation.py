"""Central simulation file.

This file serves with providing the Simulation class where the modules
are defined and set for the simulation, such as the source, the observer, 
and the propagator.

    Typical usage example:

    import rwpropa as rw

    nr_steps = 10**4
    step_size = 10**10 # [m]
    mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
    b_rms = 1 # Gaus
    b_direction = [0,0,1]
    magnetic_field = rw.OrderedBackgroundField(b_rms, b_direction).magnetic_field

    propagator = rw.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim = rw.Simulation()
    sim.add_propagator(propagator)
    sim.propagator.get_description()
    .
    .
    .
    sim.run_simulation()
    sim.save_data('data/data_sim')
"""

import numpy as np
import pandas as pd
from .source import *
from .propagator import *
from .observer import *



class Simulation():
    """Remove all step number duplicates. 

    Special observer such as the TimeEvolutionObserverLog() may introduce duplicate 
    step numbers [1,1,1,...] that lead to a wrong number of steps during the description 
    output because the step 1 will be only observed once, while it is counted several 
    times in the description output.

    Attributes:
        data: Matrix of observed data.
        source: Object of a special source type set by the user.
        observer: Object of a special observer type set by the user.
        propagator: Object of a special propagator type set by the user.
    """

    def __init__(self):
        print('start simulation')
        self.init_data()

    def init_data(self):
        self.data = [[0.0, 0.0, 0.0, 0.0, -1.0, 0.0]]
        
    def add_source(self, source):
        self.source = source
            
    def add_observer(self, observer):
        self.observer = observer

    def add_propagator(self, propagator):
        self.propagator = propagator
            
    def run_simulation(self):
        self.init_data()
        particles = self.source.particles[:]
        for p in particles:
            data_particle = p.simulate(self.observer.observer, self.propagator.propagator)
            self.data = self.data + data_particle
        self.source.reset_source()
    
    def save_data(self, file_name):
        df = pd.DataFrame(self.data[1:])
        df.columns = self.observer.get_column_names()
        df.to_pickle(file_name+".pkl")



class IsotropicSimulation():
    def __init__(self, nr_particles = 10**3, energy = 10**15, nr_steps = 10**4, diffusion_coefficient = 1.*10**21, nr_obs_steps = 600):
        self.nr_particles = nr_particles
        self.source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        self.energy = energy # [eV]
        self.nr_steps = nr_steps
        self.step_size = 1*10**12 # [m]
        self.diffusion_coefficient = diffusion_coefficient # [m^2/s]
        self.speed_of_light = 3*10**8 # [m/s]
        mfp_iso = 3*self.diffusion_coefficient/self.speed_of_light
        self.mfp = np.array([mfp_iso, mfp_iso, mfp_iso], dtype=np.float32)  # [m]
        self.nr_obs_steps = nr_obs_steps
        self.substeps = [False, False, True] # observe only steps (no substeps)
        self.sim = None

    def simulate(self, file_name = 'isotropic'):
        print('simulate with a diffusion coefficient of ', self.diffusion_coefficient, 'm²/s')
        self.sim = Simulation()
        source = PointSourceIsotropicPhi(self.energy, self.source_pos, self.nr_particles)
        self.sim.add_source(source)
        propagator = IsotropicPropagator(self.mfp, self.nr_steps, self.step_size)
        self.sim.add_propagator(propagator)
        observer = TimeEvolutionObserverLog(1, self.nr_steps, self.nr_obs_steps, self.substeps)
        self.sim.add_observer(observer)
        self.sim.run_simulation()
        self.sim.save_data(file_name)

    def get_simulation_info(self):
        if self.sim == None:
            print('run simulation first to get description of used parameters.')
        else:
            self.sim.source.get_description()
            self.sim.propagator.get_description()
            self.sim.observer.get_description()


