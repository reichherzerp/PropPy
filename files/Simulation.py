from files.Observer import Observer
import numpy as np
import pandas as pd


class Simulation():
    def __init__(self):
        print('start simulation')
        self.data = [[0.0, 0.0, 0.0, 0.0, -1.0]]
        self.add_observer(np.array([False, False, True]))
        
    def add_source(self, source):
        self.source = source
            
    def add_observer(self, observer):
        self.observer = observer
            
    def run_simulation(self, nr_steps):
        for p in self.source.particles:
            data_particle = p.simulate(self.observer, nr_steps)
            self.data = self.data + data_particle
    
    def save_data(self, file_name):
        df = pd.DataFrame(self.data[1:])
        df.columns = ['id', 'i', 'd', 'x', 'y', 'z', 'radius']
        df.to_pickle(file_name+".pkl")