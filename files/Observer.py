from numba import jit, b1, float32
import numpy as np
from numba.experimental import jitclass

simulation_spec = [
    ('observe_intermidiate', b1),
    ('observe_steps', b1[:]),
    ('pos', float32[:]),
    ('pos_prev', float32[:]),
    ('spheres', float32[:]),
]

@jitclass(simulation_spec)
class Observer():
    def __init__(self, observe_steps):
        self.observe_intermidiate = False
        self.observe_steps = observe_steps
        self.spheres = np.array([0.0], dtype=np.float32)
        print('observer')
        
    
    def observe(self, i, s, distance, pos, particle_id):
        if s == 2 and len(self.spheres) > 1:
            print('todo: implement spherical observer')
            #self.on_sphere()
        elif self.observe_steps[s]:
            if ((i == 2 or i%10 == 0) and distance > 600000000.0):
                return [particle_id, i, distance, pos[0], pos[1], pos[2], -1.0]
            else:
                return None
        else:
            return None
        
    
    def on_sphere(self, pos, pos_prev):
        #sphere = observer.on_sphere(self.pos, self.pos_prev)
        #if sphere != None:
        #print('todo')
        #particle_info.append(self.observe_particle(i))
        #particle_info.append([particle_id, i, distance, pos[0], pos[1], pos[2], sphere*1.0])
        return None