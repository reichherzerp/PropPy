from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from abc import ABC, ABCMeta, abstractmethod


observer_spec = [
    ('rms', float32),
    ('direction', float32[:]),
]

@jitclass(observer_spec)
class MagneticField():
    # base observer class that is called in the simulation by the particle class to
    # determine when to write out (observe) data. The conditions to observe can be based 
    # on the time (or step) or the coordinates of the particle.
    # - step number [unique_steps] -> time (TimeEvolutionObservers)
    # - radius of observer sphere [shperes] -> sphere around source (SphericalObservers)
    # - cartesian coordinates [box_dimensions] -> box around source (BoxObserver)
    # all special observer will create an Observer object and specify the relevant parameters
    # for the observation conditions (unique_steps, shperes, box_dimensions)

    def __init__(self, rms, direction):
        self.rms = rms
        self.direction = direction
    


class AbstractMagneticFieldMeta(ABCMeta):
    # required attributes that have to be implemented in __init__ of all
    # sub classes
    required_attributes = []

    def __call__(self, *args, **kwargs):
        # check if required attributes that have to be implemented in __init__ of all
        # sub classes are really implemented. Raise an error if not
        obj = super(AbstractMagneticFieldMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj



class AbstractMagneticField(object, metaclass=AbstractMagneticFieldMeta):
    # abstract base class for all special observers.
    # functions with the label @abstractmethod have to be implemented in 
    # the special observer classes

    # all required_attributes have to be implemented in sub classes
    required_attributes = ['rms', 'direction']
 
    @abstractmethod
    def __init__(self, order):
        # implementation required in all sub classes.
        # all required_attributes have to be implemented in sub classes
        pass



class OrderedBackgroundField(AbstractMagneticField):

    def __init__(self, rms, direction):
        self.rms = rms
        self.direction = direction
        self.magnetic_field = MagneticField(rms, np.array(direction, dtype=np.float32))
        #self.init_observer(substeps)





