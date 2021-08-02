from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass

simulation_spec = [
    ('isotropic', b1),
    ('nr_steps', int32),
    ('step_size', float32)
]

@jitclass(simulation_spec)
class Propagator():
    def __init__(self, nr_steps, step_size):
        print('Propagator initialized')
        self.isotropic = False
        self.nr_steps = nr_steps
        self.step_size = step_size