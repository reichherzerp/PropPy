from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass

simulation_spec = [
    ('observe_intermidiate', b1),
    ('substeps', b1[:]),
    ('all_steps', b1),
    ('pos', float32[:]),
    ('steps', int32[:]),
    ('pos_prev', float32[:]),
    ('spheres', float32[:]),
]

@jitclass(simulation_spec)
class Observer():
    def __init__(self, steps, substeps):
        self.observe_intermidiate = False
        self.substeps = substeps
        self.spheres = np.array([0.0], dtype=np.float32)
        self.all_steps = False
        if -1 in steps:
            # the -1 is the key to say that all steps should be observed
            self.all_steps = True
        self.steps = steps
        print('observer')
        
    
    def observe(self, i, substep, distance, pos, particle_id):
        if substep == 2 and len(self.spheres) > 1:
            print('todo: implement spherical observer')
            #self.on_sphere()
        elif self.substeps[substep]:
            if self.all_steps or i in self.steps:
                return [particle_id, i, distance, pos[0], pos[1], pos[2], -1.0, substep]
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

class TimeEvolutionObserverLog():
    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        substeps_bool = np.array(substeps) 
        steps = np.logspace(np.log10(min_steps), np.log10(max_steps), nr_steps)
        steps_int32 = np.array(steps, dtype=np.int32) # np.array([-1])

        self.observer = Observer(steps_int32, substeps_bool)

class TimeEvolutionObserverLin():
    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        substeps_bool = np.array(substeps) 
        steps = np.linspace(min_steps, max_steps, nr_steps)
        steps_int32 = np.array(steps, dtype=np.int32) # np.array([-1])

        self.observer = Observer(steps_int32, substeps_bool)

class TimeEvolutionObserver():
    def __init__(self, steps, substeps):
        substeps_bool = np.array(substeps) 
        steps_int32 = np.array(steps, dtype=np.int32) # np.array([-1])

        self.observer = Observer(steps_int32, substeps_bool)