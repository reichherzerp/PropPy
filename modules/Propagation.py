import random
from numba import jit, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass
import numpy as np

spec = [   
    ('step_size', float32), # a simple scalar field
]

@jitclass(spec)
class Propagation():
    def __init__(self, step_size):
        self.step_size = step_size
        
    def move(self, particle):
        kappa = particle.diffusion_tensor
        tau_step = 1.0
        chi = np.array([1.0, 1.0, 1.0])
        v = 1.0
        chi = chi/np.sum(chi)*v
        tau = kappa*tau_step**2/chi**2
        xi = 2/tau
        prop_turn = tau_step*xi
        random_value = random.random()
        
        step_size = 1
        p = 1.0*step_size
        particle.t = particle.t + 1*step_size
        
        
        ### 1. change direction (if needed)
        if particle.t > particle.gyro_radius:
            # after gyroradius ~ half gyroorbit, particle 
            # needs to change direction perp to background field
            if (random_value <= prop_turn[0]):
                particle.direction[0] = particle.direction[0]*(-1)
            if (random_value <= prop_turn[1]):
                particle.direction[1] = particle.direction[1]*(-1)
        if (random_value <= prop_turn[2]):
            particle.direction[2] = particle.direction[2]*(-1)
          
        ### 2. move along direction
        abs_direction = (particle.direction[0]**2+particle.direction[1]**2+particle.direction[2]**2)**0.5
        speed = particle.direction / abs_direction * v
        d = speed * tau_step 
        particle.pos[0] = particle.pos[0] + d[0]
        particle.pos[1] = particle.pos[1] + d[1]
        particle.pos[2] = particle.pos[2] + d[2]
        return particle