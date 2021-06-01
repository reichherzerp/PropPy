from numba.experimental import jitclass
import matplotlib.pyplot as plt
from numba import jit, int32, float32

source_spec = [
    ('nr_particles', int32), # a simple scalar field
    ('position', float32[:]), # an array field  
    ('diffusion_tensor', float32[:]),
    ('gyro_radius', float32),
]

@jitclass(source_spec)
class Source():
    def __init__(self, nr_particles, position, gyro_radius, diffusion_tensor):
        self.nr_particles = nr_particles
        self.position = position
        self.gyro_radius = gyro_radius
        self.diffusion_tensor = diffusion_tensor