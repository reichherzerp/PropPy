import random
import numpy as np
from numba import jit, int32, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass
import numpy as np
from modules.Particle import Particle
from modules.Propagation import Propagation
from modules.Observer import Observer


simulation_spec = [
    ('particles', types.ListType(Particle.class_type.instance_type)),
    ('propagation', Propagation.class_type.instance_type),
    ('observer', Observer.class_type.instance_type),
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

    def add_observer(self, observer):
        self.observer = observer
    
    def run_simulation(self):
        id = []
        x = []
        y = []
        z = []
        time = []
        radius = []
        
        for i, t in enumerate(self.time):
            for p in self.particles:

                self.propagation.move_polar(p)
                
                ### observer
                observe = self.observer.observe(i, p)
                should_observe = observe[0]

                ### output data
                if should_observe:
                    id.append(p.id)
                    x.append(p.pos[0])
                    y.append(p.pos[1])
                    z.append(p.pos[2])
                    time.append(t)
                    radius.append(observe[1])

        return id, time, x, y, z, radius