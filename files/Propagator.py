from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass

simulation_spec = [
    ('isotropic', b1)
]

@jitclass(simulation_spec)
class Propagator():
    def __init__(self):
        print('Propagator initialized')
        self.isotropic = False