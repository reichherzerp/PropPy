import numpy as np
from numba import jit, float32, int32, types, typed
from numba.typed import List
from numba.experimental import jitclass



spec = [
    ('diffusion_tensor', float32[:]), 
    ('gyro_radius', float32), # a simple scalar field 
    ('t', float32),
    ('pos', float32[:]),          # an array field
    ('pos_previous', float32[:]),
    ('direction', int32[:]),
    ('id', int32),
    ('phi', float32), 
]

@jitclass(spec)
class Particle():
    def __init__(self, id, gyro_radius, diffusion_tensor):
        self.id = id
        self.gyro_radius = gyro_radius
        self.diffusion_tensor = diffusion_tensor
        self.t = 0
        self.phi = 0
        self.pos = np.zeros(3, dtype=np.float32)
        self.pos_previous = np.zeros(3, dtype=np.float32)
        
    def setDirection(self, direction):
        self.direction = direction
        
    def kappa(self, axis):
        return self.pos[axis]**2/(2.0)