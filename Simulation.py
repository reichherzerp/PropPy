import matplotlib.pyplot as plt
import random
import numpy as np
from numba import jit, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass
import numpy as np
from .Particle import Particle

simulation_spec = [
    ('particles', types.ListType(Particle.class_type.instance_type)),
    ('propagation', Propagation.class_type.instance_type),
    ('time', float32[:]),
]

@jitclass(simulation_spec)
class Simulation():
    def __init__(self):
        print('start simulation')
  
    def addParticles(self, source):
        particles = List()
        for j in range(source.nr_particles):
            particle = Particle(source.gyro_radius, source.free_mean_path_para, source.free_mean_path_perp)
            direction = List()
            random_values = np.random.randint(low=0, high=1, size=3)*2-1
            [direction.append(x) for x in random_values]
            particle.setDirection(direction)
            particles.append(particle)
            
        self.particles = particles

    def addPropagation(self, propagation, time):
        self.propagation = propagation
        self.time = time
        
    def distribution(self, axis):
        data = []
        for p in self.particles:
            data.append(p.pos[axis])
        return data
    
    def runSimulation(self, time):
        kappa_perp = []
        kappa_para = []
        x = [0]
        y = [0]
        for i, t in enumerate(time):
            kappa_perp_sum = 0
            kappa_para_sum = 0
            particles = List()
            [particles.append(self.propagation.move(p)) for p in self.particles]
            self.particles = particles