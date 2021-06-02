from numba import jit, int32, float32, types, typed
from numba.typed import List
from numba import numba
from numba.experimental import jitclass


observer_spec = [
    ('observer_type', int32),
    ('sphere_radii', float32[:]),
]

@jitclass(observer_spec)
class Observer():
    def __init__(self, observer_type, sphere_radii):
        self.observer_type = observer_type
        self.sphere_radii = sphere_radii 
        print('init observer completed')