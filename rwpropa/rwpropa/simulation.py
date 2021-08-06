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
        self.observer = observer

    def add_propagator(self, propagator):
        rms = 10
        direction = [0,0,1]
        #self.magnetic_field = OrderedBackgroundField(rms, direction)
        #propagator.set_magnetic_field(self.magnetic_field.magnetic_field)
        self.propagator = propagator.propagator
            
    def run_simulation(self):
        self.init_data()
        particles = self.source.particles[:]
        for p in particles:
            data_particle = p.simulate(self.observer.observer, self.propagator)
            self.data = self.data + data_particle
        self.source.reset_source()
    
    def save_data(self, file_name):
        df = pd.DataFrame(self.data[1:])
        df.columns = self.observer.get_column_names()
        df.to_pickle(file_name+".pkl")