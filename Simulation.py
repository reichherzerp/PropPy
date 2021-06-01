import random
import numpy as np
from numba import jit, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass
import numpy as np
from modules.Particle import Particle
from modules.Propagation import Propagation
from modules.Source import Source

simulation_spec = [
    ('particles', types.ListType(Particle.class_type.instance_type)),
    ('propagation', Propagation.class_type.instance_type),
    ('time', float32[:]),
]

@jitclass(simulation_spec)
class Simulation():
    def __init__(self):
        print('init simulation completed')
  
    def addParticles(self, source):
        particles = List()
        for j in range(source.nr_particles):
            particle = Particle(source.gyro_radius, source.diffusion_tensor)
            #direction = List()
            random_values = np.ones(3, dtype=np.int32)
            random_values[0] = random.randint(0, 1)*2-1
            random_values[1] = random.randint(0, 1)*2-1
            random_values[2] = random.randint(0, 1)*2-1
            particle.setDirection(random_values)
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
    
    def runSimulation(self):
        kappa_perp = []
        kappa_para = []
        x = [0]
        y = [0]
        for i, t in enumerate(self.time):
            kappa_perp_sum = 0
            kappa_para_sum = 0
            particles = List()
            [particles.append(self.propagation.move(p)) for p in self.particles]
            for p in particles:
                kappa_perp_sum = kappa_perp_sum + (p.kappa(0)+p.kappa(1))/2
                kappa_para_sum = kappa_para_sum + p.kappa(2)
    
            kappa_para.append(kappa_para_sum/len(particles)/t)
            kappa_perp.append(kappa_perp_sum/len(particles)/t)
            self.particles = particles
        return [kappa_para, kappa_perp]