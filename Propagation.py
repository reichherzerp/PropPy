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
        step_size = 1
        p = 1.0*step_size
        particle.t = particle.t + 1*step_size
        
        if particle.t > particle.gyro_radius:
            # after gyroradius ~ half gyroorbit, particle 
            # needs to change direction perp to background field
            if (random.randint(0, particle.free_mean_path_perp) == 1):
                particle.direction[0] = particle.direction[0]*(-1)
            if (random.randint(0, particle.free_mean_path_perp) == 0):
                particle.direction[1] = particle.direction[1]*(-1)
        if particle.t > particle.free_mean_path_para:
            # after mean-free length, particle can also change 
            # direction parallel to mean field
            if (random.randint(0, particle.free_mean_path_para) == 0):
                particle.direction[2] = particle.direction[2]*(-1)
            
        absDirection = (particle.direction[0]**2+particle.direction[1]**2+particle.direction[2]**2)**0.5
        particle.pos[0] = particle.pos[0] + p*particle.direction[0]/absDirection
        particle.pos[1] = particle.pos[1] + p*particle.direction[1]/absDirection
        particle.pos[2] = particle.pos[2] + p*particle.direction[2]/absDirection
        return particle