import random
import numpy as np
from numba import jit, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass
import numpy as np
from modules.Particle import Particle
from modules.Propagation import Propagation


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
            particle = Particle(j, source.gyro_radius, source.diffusion_tensor)
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

    def addOutput(self, output):
        self.output = output
        
    def distribution(self, axis):
        data = []
        for p in self.particles:
            data.append(p.pos[axis])
        return data
    
    def runSimulation(self):
        id = []
        x = []
        y = []
        z = []
        time = []
        
        for i, t in enumerate(self.time):
            particles = List()
            for p in self.particles:
                particles.append(self.propagation.move(p))
                ### observer
                if i < 1000 or i % 1000 == 0:
                    id.append(p.id)
                    x.append(p.pos[0])
                    y.append(p.pos[1])
                    z.append(p.pos[2])
                    time.append(t)
   
            self.particles = particles
  
        return id, time, x, y, z

        
