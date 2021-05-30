import numpy as np
from numba import jit, float32, types, typed
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
class Particle():
    def __init__(self, gyro_radius, free_mean_path_para, free_mean_path_perp):
        self.gyro_radius = gyro_radius
        self.free_mean_path_para = free_mean_path_para
        self.free_mean_path_perp = free_mean_path_perp
        self.t = 0
        self.pos = np.zeros(3, dtype=np.float32)
        
    def setDirection(self, direction):
        self.direction = direction
        
    def kappa(self, axis):
        return self.pos[axis]**2/(2*self.t)