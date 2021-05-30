import matplotlib.pyplot as plt
import random
import numpy as np
from numba import jit, int32, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass

spec = [
    ('free_mean_path_para', float32), # a simple scalar field
    ('free_mean_path_perp', float32), 
    ('gyro_radius', float32), 
    ('t', float32),
    ('pos', float32[:]),          # an array field
    ('direction', types.ListType(types.int64))     
]

@jitclass(spec)
class Popagation():
    def __init__(self, gyro_radius, free_mean_path_para, free_mean_path_perp, direction):
        self.gyro_radius = gyro_radius
        self.free_mean_path_para = free_mean_path_para
        self.free_mean_path_perp = free_mean_path_perp
        self.t = 0
        self.pos = np.zeros(3, dtype=np.float32)
        self.direction = direction
                 
     
    def move(self, step_size):
        p = 1.0*step_size
        self.t = self.t + 1*step_size
        
        if self.t > self.gyro_radius:
            # after gyroradius ~ half gyroorbit, particle 
            # needs to change direction perp to background field
            if (random.randint(0, self.free_mean_path_perp) == 1):
                self.direction[0] = self.direction[0]*(-1)
            if (random.randint(0, self.free_mean_path_perp) == 0):
                self.direction[1] = self.direction[1]*(-1)
        if self.t > self.free_mean_path_para:
            # after mean-free length, particle can also change 
            # direction parallel to mean field
            if (random.randint(0, self.free_mean_path_para) == 0):
                self.direction[2] = self.direction[2]*(-1)
            
        absDirection = (self.direction[0]**2+self.direction[1]**2+self.direction[2]**2)**0.5
        self.pos[0] = self.pos[0] + p*self.direction[0]/absDirection
        self.pos[1] = self.pos[1] + p*self.direction[1]/absDirection
        self.pos[2] = self.pos[2] + p*self.direction[2]/absDirection
        
    def kappa(self, axis):
        return self.pos[axis]**2/(2*self.t)