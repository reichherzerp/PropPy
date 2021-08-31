"""Customized ordered magnetic fields.

Ordered magnetic fields determine the directions of the parallel and perpendicular 
diffusion coefficients of the diffusion tensor. The parallel diffusion coefficient
is used for the parallel direction to the ordered background magnetic field. 
The magnetic fields are used by the propagation class to determine the orientation 
of the field and how to align the diffusion tensor. 

    Typical usage example:

    import rwpropa as rw

    nr_steps = 10**4
    step_size = 10**10 # [m]
    mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
    b_rms = 1 # Gaus
    b_direction = [0,0,1]
    magnetic_field = rw.OrderedBackgroundField(b_rms, b_direction).magnetic_field

    propagator = rw.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim = rw.Simulation()
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
    """ Base observer class that is called in the simulation by the particle class to
    determine when to write out (observe) data. 
     
    The conditions to observe can be based 
    on the time (or step) or the coordinates of the particle.
    - step number [unique_steps] -> time (TimeEvolutionObservers)
    - radius of observer sphere [shperes] -> sphere around source (SphericalObservers)
    - cartesian coordinates [box_dimensions] -> box around source (BoxObserver)
    all special observer will create an Observer object and specify the relevant parameters
    for the observation conditions (unique_steps, shperes, box_dimensions)

    Attributes:
        rms: A float32 indicating the root-mean square value of the magnetic field.
        direction: An float32 array indicating the direction of the ordered magnetic field.
    """

    def __init__(self, rms, direction):
        self.rms = rms
        self.direction = direction
    


class AbstractMagneticFieldMeta(ABCMeta):
    """ Abstract meta class to check if all required attributes are implemented in the 
    sub classes.
    """
    required_attributes = []

    def __call__(self, *args, **kwargs):
        """ Checks if required attributes that have to be implemented in __init__ of all
        sub classes are really implemented. 

        Raises:
            ValueError: an error if not all required attributes are implemented.
        """
        obj = super(AbstractMagneticFieldMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj



class AbstractMagneticField(object, metaclass=AbstractMagneticFieldMeta):
    """Abstract base class for all special observers.
    
    Functions with the label @abstractmethod have to be implemented in the special 
    observer classes.

    Attributes:
        rms: A float32 indicating the root-mean square value of the magnetic field.
        direction: An float32 array indicating the direction of the ordered magnetic field.
        required_attributes: A list with attributes that have to be implemented in sub classes.
    """

    required_attributes = [
        'rms', 
        'direction'
    ]
 
    @abstractmethod
    def __init__(self, order):
        """ Implementation required in all sub classes. All required_attributes have to 
        be implemented in sub classes.
        """
        pass



class OrderedBackgroundField(AbstractMagneticField):
    """Ordered background magnetic field.

    The user can specify the root-mean square value and the direction of the ordered 
    magnetic field. The magnetic field is static and points everywhere in the specifyed
    direction with the rms vlaue that is specifyed.

    Attributes:
        rms: A float32 indicating the root-mean square value of the magnetic field.
        direction: An float32 array indicating the direction of the ordered magnetic field.
        magnetic_field: An object of the @jitclass MagneticField that can be passed to other jitclasses.
    """

    def __init__(self, rms, direction):
        """Initializes the required parameters and creates the @jitclass magnetic field
        with the user specifyed rms and direction."""
        self.rms = rms
        self.direction = np.array(direction, dtype=np.float32)
        self.magnetic_field = MagneticField(self.rms, self.direction)



class DefaultBackgroundField(AbstractMagneticField):
    """Ordered default background magnetic field.

    The user can only specify the root-mean square value of the ordered  magnetic field. 
    The magnetic field is static and points everywhere along the z-axis with the rms vlaue 
    that is specifyed.

    Attributes:
        rms: A float32 indicating the root-mean square value of the magnetic field.
        magnetic_field: An object of the @jitclass MagneticField that can be passed to other jitclasses.
    """

    def __init__(self, rms):
        """Initializes the required parameters and creates the @jitclass magnetic field
        with the user specifyed rms and the default direction along the z-axis."""
        self.rms = rms
        self.direction = np.array([0,0,1], dtype=np.float32)
        self.magnetic_field = MagneticField(self.rms, self.direction)