from numba import jit, int32, float32, types, typed
from numba.typed import List
import numpy as np
from numba.experimental import jitclass


observer_spec = [
    ('observer_type', int32),
    ('should_observe_now', int32),
    ('sphere_radii', float32[:]),
]

@jitclass(observer_spec)
class Observer():
    def __init__(self):
        self.observer_type = 0
        self.should_observe_now = 0
        self.sphere_radii = np.array([-1.0], dtype=np.float32)
        print('init observer completed')

    def add_observer_spheres(self, sphere_radii):
        self.sphere_radii = sphere_radii 
        self.observer_type = 1

    def observe(self, i, p):
        if self.observer_type == 0:
            if i < 1000 or i % 1000 == 0:
                return True, -1.0
        if self.observer_type == 1:
            r2 = p.pos[0]**2+p.pos[1]**2+p.pos[2]**2
            r2_previous = p.pos_previous[0]**2+p.pos_previous[1]**2+p.pos_previous[2]**2
            for r_sphere in self.sphere_radii:
                r2_sphere = r_sphere**2
                if (r2 >= r2_sphere and r2_previous < r2_sphere) or (r2 <= r2_sphere and r2_previous > r2_sphere):
                    return True, r_sphere 
                    ### in this case, the particles crossed the sphere with radius^2 = r2_sphere
                
        return False, -1.0

        