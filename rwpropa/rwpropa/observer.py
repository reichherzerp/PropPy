from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from abc import ABC, ABCMeta, abstractmethod
from .particle_state import *


observer_spec = [
    ('substeps', b1[:]),
    ('all_steps', b1),
    ('pos', float32[:]),
    ('steps', int32[:]),
    ('pos_prev', float32[:]),
    ('spheres', float32[:]),
    ('box_dimensions', float32[:]),
    ('ps', ParticleState.class_type.instance_type),
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
        self.spheres = np.array([-1.0, 10**5], dtype=np.float32)
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
        

    def observe(self, ps):
        # decide if the current particle state should be observed based on the criterions specified 
        # in the observer instance
        if ps.substep == 2 and len(self.spheres) > 1:
            on_sphere = self.check_on_sphere(ps)
            if on_sphere != None:
                return on_sphere
        if self.substeps[ps.substep]:
            if self.all_steps or ps.step in self.steps:
                return self.data_row(ps, -1.0)
            else:
                return None
        else:
            return None


    def check_on_sphere(self, ps):
        radius = 0
        for i in range(ps.dimensions):
            radius = radius + ps.pos[i]**2

        if radius**0.5 > self.spheres[1]:
            return self.data_row(ps, self.spheres[1])
        else: 
            return None


    def data_row(self, ps, radius):
        radius = radius
        data_row_list = [
            ps.particle_id, 
            ps.step, 
            ps.distance, 
            ps.pos[0], 
            ps.pos[1], 
            ps.pos[2], 
            ps.phi, 
            ps.pitch_angle, 
            radius, 
            ps.substep
        ]
        return data_row_list
 

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



#-------------------------------------------------------------------------------
# below are the abstract base class and all sub classes of the special observers
# that have to be added to the simulation. Each special observer stores a
# Observer object in its instance parameter to be used in the simulation.
# This diversions is needed because numba does not support 
# inheritance via ABC and Propagator() needs the label @jitclass as it is called 
# during the numba optimized simulation loop of the run_simulation() function. 
# This workaround supports both concepts with the 
# advantages of fast code and easy addition of new observers where the structure 
# is now defined by the Abstract Base class and enforeced via the ABCMeta class



class AbstractSpecialObserverMeta(ABCMeta):
    # required attributes that have to be implemented in __init__ of all
    # sub classes
    required_attributes = []

    def __call__(self, *args, **kwargs):
        # check if required attributes that have to be implemented in __init__ of all
        # sub classes are really implemented. Raise an error if not
        obj = super(AbstractSpecialObserverMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj



class AbstractSpecialObserver(object, metaclass=AbstractSpecialObserverMeta):
    # abstract base class for all special observers.
    # functions with the label @abstractmethod have to be implemented in 
    # the special observer classes

    # all required_attributes have to be implemented in sub classes
    required_attributes = [
        'steps', 
        'substeps_bool'
    ]
 
    @abstractmethod
    def __init__(self, order):
        # implementation required in all sub classes.
        # all required_attributes have to be implemented in sub classes
        pass

    def init_observer(self, substeps):
        # set the important parameters and call the @abstractmethods that are implemented
        # in each special observer class that are derived from the current abstract base class
        self.column = ['id', 'i', 'd', 'x', 'y', 'z', 'phi', 'pitch_angle', 'radius', 'sub_step']
        self.substeps_bool = np.array(substeps)
        self.steps = self.set_steps_int() 
        # have to store all relevant observation parameters in the Observer class that 
        # has the @jitclass label from numba. This is important, as the Particle class is also 
        # labeled with @jitclass and can thus only call @jitclass classes. This is for
        # performance resaons.
        self.observer = Observer(self.steps, self.substeps_bool)

    @abstractmethod
    def set_steps(self):
        # set the number of steps that should be observed
        pass

    def set_steps_int(self):
        # convert list of steps to np array of ints that is known in the @jitclass Observer()
        steps = self.set_steps()
        steps_int32 = np.array(steps, dtype=np.int32)
        return steps_int32
    
    @abstractmethod
    def get_description_observer_type(self):
        # give the name of each special observer
        pass

    def get_description(self):
        self.observer.get_description_general()
        self.get_description_observer_type()
        self.observer.get_description_parameters()

    def get_column_names(self):
        return self.column

    

class ObserverAllSteps(AbstractSpecialObserver):

    def __init__(self, substeps):
        self.init_observer(substeps)

    def set_steps(self):
        steps = [-1]
        return steps

    def get_description_observer_type(self):
        print('observer tpye: ObserverAllSteps')



class TimeEvolutionObserverLog(AbstractSpecialObserver):

    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        self.min_steps = min_steps
        self.max_steps = max_steps
        self.nr_steps = nr_steps

        self.init_observer(substeps)
        
    def set_steps(self):
        steps = np.logspace(np.log10(self.min_steps), np.log10(self.max_steps), self.nr_steps)
        return steps

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserverLog')



class TimeEvolutionObserverLin(AbstractSpecialObserver):

    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        self.min_steps = min_steps
        self.max_steps = max_steps
        self.nr_steps = nr_steps

        self.init_observer(substeps)

    def set_steps(self):
        steps = np.linspace(self.min_steps, self.max_steps, self.nr_steps)
        return steps

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserverLin')
    


class TimeEvolutionObserver(AbstractSpecialObserver):

    def __init__(self, steps_input, substeps):
        self.steps_input = steps_input

        self.init_observer(substeps)

    def set_steps(self):
        steps = self.steps_input 
        return steps

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserver')