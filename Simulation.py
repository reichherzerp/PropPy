import random
import numpy as np
from numba import jit, int32, float32, types, typed
from numba.typed import List
import numba
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
  
    def add_particles(self, source):
        particles = List()
        for j in range(source.nr_particles):
            particle = Particle(j, source.gyro_radius, source.diffusion_tensor)
            random_values = np.ones(3, dtype=np.int32)
            random_values[0] = random.randint(0, 1)*2-1
            random_values[1] = random.randint(0, 1)*2-1
            random_values[2] = random.randint(0, 1)*2-1
            particle.setDirection(random_values)
            particles.append(particle)
            
        self.particles = particles

    def add_propagation(self, propagation, time):
        self.propagation = propagation
        self.time = time
        
    def distribution(self, axis):
        data = []
        for p in self.particles:
            data.append(p.pos[axis])
        return data
    
    def run_simulation(self, observer_type):
        id = []
        x = []
        y = []
        z = []
        time = []
        
        for i, t in enumerate(self.time):
            for p in self.particles:
                self.propagation.move(p)
                ### observer
                if observer_type == 0:
                    if i < 1000 or i % 1000 == 0:
                        id.append(p.id)
                        x.append(p.pos[0])
                        y.append(p.pos[1])
                        z.append(p.pos[2])
                        time.append(t)
                if observer_type == 1:
                    r2_sphere = (1*10**2)**2
                    r2 = p.pos[0]**2+p.pos[1]**2+p.pos[2]**2
                    r2_previous = p.pos_previous[0]**2+p.pos_previous[1]**2+p.pos_previous[2]**2
                    #print(r2, r2_previous, r2_sphere)
                    if (r2 >= r2_sphere and r2_previous < r2_sphere) or (r2 <= r2_sphere and r2_previous > r2_sphere):
                        ### in this case, the particles crossed the sphere with radius^2 = r2_sphere
                        id.append(p.id)
                        x.append(p.pos[0])
                        y.append(p.pos[1])
                        z.append(p.pos[2])
                        time.append(t)

        return id, time, x, y, z