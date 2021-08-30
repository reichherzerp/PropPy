"""Customized ordered magnetic fields.

Ordered magnetic fields determine the directions of the parallel and perpendicular 
diffusion coefficients of the diffusion tensor. The parallel diffusion coefficient
is used for the parallel direction to the ordered background magnetic field. 
The magnetic fields are used by the propagation class to determine the orientation 
of the field and how to align the diffusion tensor. 

    Typical usage example:

    nr_steps = 10**4
    step_size = 10**10 # [m]
    mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
    b_rms = 1 # Gaus
    b_direction = [0,0,1]
    magnetic_field = rw.OrderedBackgroundField(b_rms, b_direction).magnetic_field

    propagator = rw.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim.add_propagator(propagator)
    sim.propagator.get_description()

"""

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
    required_attributes = [
        'rms', 
        'direction'
    ]
 
    @abstractmethod
    def __init__(self, order):
        # implementation required in all sub classes.
        # all required_attributes have to be implemented in sub classes
        pass



class OrderedBackgroundField(AbstractMagneticField):

    def __init__(self, rms, direction):
        self.rms = rms
        self.direction = np.array(direction, dtype=np.float32)
        self.magnetic_field = MagneticField(self.rms, self.direction)



class DefaultBackgroundField(AbstractMagneticField):

    def __init__(self, rms):
        self.rms = rms
        self.direction = np.array([0,0,1], dtype=np.float32)
        self.magnetic_field = MagneticField(self.rms, self.direction)





