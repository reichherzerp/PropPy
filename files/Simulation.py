import pandas as pd


class Simulation():
    def __init__(self):
        print('start simulation')
        self.init_data()

    def init_data(self):
        self.data = [[0.0, 0.0, 0.0, 0.0, -1.0, 0.0]]
        
    def add_source(self, source):
        self.source = source
            
    def add_observer(self, observer):
        self.observer = observer.observer

    def add_propagator(self, propagator):
        self.propagator = propagator
            
    def run_simulation(self):
        if len(self.data) > 1:
            self.init_data()
        print('init source')
        self.source.init_source()
        particles = self.source.particles[:]
        for p in particles:
            data_particle = p.simulate(self.observer, self.propagator)
            self.data = self.data + data_particle
    
    def save_data(self, file_name):
        df = pd.DataFrame(self.data[1:])
        df.columns = ['id', 'i', 'd', 'x', 'y', 'z', 'radius', 'step', 'phi', 'pitch_angle']
        df.to_pickle(file_name+".pkl")