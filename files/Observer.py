from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from abc import ABC, abstractmethod


observer_spec = [
    ('substeps', b1[:]),
    ('all_steps', b1),
    ('pos', float32[:]),
    ('steps', int32[:]),
    ('pos_prev', float32[:]),
    ('spheres', float32[:]),
    ('box_dimensions', float32[:]),
]

@jitclass(observer_spec)
class Observer():
    # base observer class that is called in the simulation by the particle class to
    # determine when to write out (observe) data. The conditions to observe can be based 
    # on the time (or step) or the coordinates of the particle.
    # - step number [unique_steps] -> time (TimeEvolutionObservers)
    # - radius of observer sphere [shperes] -> sphere around source (SphericalObservers)
    # - cartesian coordinates [box_dimensions] -> box around source (BoxObserver)
    # all special observer will create an Observer object and specify the relevant parameters
    # for the observation conditions (unique_steps, shperes, box_dimensions)
    def __init__(self, steps, substeps):
        self.substeps = substeps
        self.spheres = np.array([0.0], dtype=np.float32)
        # distance of the box boundaries in cartesian coords. [x, y, z]
        self.box_dimensions = np.array([0.0], dtype=np.float32)
        self.all_steps = False
        if -1 in steps:
            # the -1 is the key to say that all steps should be observed
            self.all_steps = True
        self.steps = self.get_unique_steps(steps)
        print('number steps: ', len(self.steps))
        print('Observer initialized')
    

    def get_unique_steps(self, steps):
        # remove all step number duplicates. Special observer such as the 
        # TimeEvolutionObserverLog() may introduce duplicate step numbers [1,1,1,...] 
        # that lead to a wrong number of steps during the description output because the step 1
        # will be only observed once, while it is counted several times in the description output
        unique_steps = []
        for s in steps:
            if s not in unique_steps:
                unique_steps.append(s)
        return np.array(unique_steps, dtype=np.int32)
        

    def observe(self, i, substep, distance, pos, particle_id, phi, pitch_angle):
        if substep == 2 and len(self.spheres) > 1:
            print('todo: implement spherical observer')
            #self.on_sphere()
        elif self.substeps[substep]:
            if self.all_steps or i in self.steps:
                radius = -1.0 # default
                return [particle_id, i, distance, pos[0], pos[1], pos[2], phi, pitch_angle, radius, substep]
            else:
                return None
        else:
            return None
 

    def get_description(self):
        # note: description does not contain the information of the underling special observer 
        # (if there was one)
        # that was used during the observer initialization. To get this info, get_description(self) 
        # has to be called directly on the instance of the special observer class (see tutorials for details)
        self.get_description_general()
        self.get_description_parameters()


    def get_description_general(self):
        # called by all special observer classes below.
        # introduction of the description output
        print("""Description Observer:
                The observer defines the conditions for when to write data to the output.\n""")


    def get_description_parameters(self):   
        # called by all special observer classes below.
        # print out all relevant instance parameters
        print('steps [0:10]: ' , self.steps[0:10])
        print('steps [-11:-1]: ' ,self.steps[-11:-1])
        print('nr steps: ' , len(self.steps))
        print('substeps: ', self.substeps)  
        print('all_steps: ', self.all_steps) 



class ObserverData():
    def __init__(self):
        self.column = ['id', 'i', 'd', 'x', 'y', 'z', 'phi', 'pitch_angle', 'radius', 'sub_step']


class AbstractSpecialObserver(ABC):
 
    def __init__(self):
        #self.value = value
        super().__init__()

    @abstractmethod
    def set_steps(self):
        # set the number of steps that should be observed
        pass
    
    @abstractmethod
    def get_description_observer_type(self):
        # give the name of each special observer
        pass

    def get_description(self):
        self.observer.get_description_general()
        self.get_description_observer_type()
        self.observer.get_description_parameters()

    def get_column_names(self):
        return ObserverData().column

    


class ObserverAllSteps(AbstractSpecialObserver):
    def __init__(self, substeps):
        self.set_steps()
        substeps_bool = np.array(substeps) 
        self.observer = Observer(self.steps, substeps_bool)

    def set_steps(self):
        steps = [-1]
        steps_int32 = np.array(steps, dtype=np.int32)
        self.steps = steps_int32

    def get_description_observer_type(self):
        print('observer tpye: ObserverAllSteps')



class TimeEvolutionObserverLog(AbstractSpecialObserver):
    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        self.min_steps = min_steps
        self.max_steps = max_steps
        self.nr_steps = nr_steps
        substeps_bool = np.array(substeps) 
        self.set_steps()
        self.observer = Observer(self.steps, substeps_bool)

    def set_steps(self):
        steps = np.logspace(np.log10(self.min_steps), np.log10(self.max_steps), self.nr_steps)
        steps_int32 = np.array(steps, dtype=np.int32)
        self.steps = steps_int32

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserverLog')



class TimeEvolutionObserverLin(AbstractSpecialObserver):
    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        self.min_steps = min_steps
        self.max_steps = max_steps
        self.nr_steps = nr_steps
        self.set_steps()
        substeps_bool = np.array(substeps) 
        self.observer = Observer(self.steps, substeps_bool)

    def set_steps(self):
        steps = np.linspace(self.min_steps, self.max_steps, self.nr_steps)
        steps_int32 = np.array(steps, dtype=np.int32) 
        self.steps = steps_int32

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserverLin')
    


class TimeEvolutionObserver(AbstractSpecialObserver):
    def __init__(self, steps, substeps):
        self.steps = steps
        substeps_bool = np.array(substeps) 
        self.set_steps()
        self.observer = Observer(self.steps, substeps_bool)

    def set_steps(self):
        steps = np.linspace(self.min_steps, self.max_steps, self.nr_steps)
        steps_int32 = np.array(steps, dtype=np.int32) 
        self.steps = steps_int32

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserver')