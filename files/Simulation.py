from files.Observer import Observer
from files.Particle import Particle
import numpy as np
import pandas as pd

class Simulation():
    def __init__(self):
        print('start simulation')
        self.data = [[0.0, 0.0, 0.0, 0.0, -1.0]]
        self.particles = []
        self.add_observer(np.array([False, False, True]))
        
    def add_source(self, gyro_radius, mean_free_path, pos, nr_particles):
        for i in range(nr_particles):
            particle_id = i
            p = Particle(particle_id, gyro_radius, mean_free_path, pos)
            self.particles.append(p)
            
    def add_observer(self, observe_substeps):
        self.observer = Observer(observe_substeps)
            
    def run_simulation(self, nr_steps):
        for p in self.particles:
            data_particle = p.simulate(self.observer, nr_steps)
            self.data = self.data + data_particle
    
    def save_data(self, file_name):
        df = pd.DataFrame(self.data[1:])
        df.columns = ['id', 'i', 'd', 'x', 'y', 'z', 'radius']
        df.to_pickle(file_name+".pkl")