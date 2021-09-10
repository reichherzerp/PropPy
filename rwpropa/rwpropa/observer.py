"""The Observer determines during the simulation when and what data to write out (observe).

In each simulation step, the current particle state is evaluated by the Observer to check
if one of the observing conditions is satisfied. The conditions to observe can be based 
on the time (-> step) or the coordinates of the particle.

    Typical usage example:

    import rwpropa as rw
    
    substeps = [False, False, True] # observe only steps (no substeps)
    spheres = [1*10**13, 5*10**12]

    observer = rw.SphericalObserver(substeps, spheres)

    sim = rw.Simulation()
    sim.add_observer(observer)
    sim.observer.get_description()
"""


from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from abc import ABC, ABCMeta, abstractmethod
from .particle_state import *


observer_spec = [
    ('substeps', b1[:]),
    ('all_steps', b1),
    ('steps', int32[:]),
    ('pos', float32[:]),
    ('spheres', float32[:]),
    ('ps', ParticleState.class_type.instance_type),
]

@jitclass(observer_spec)
class Observer():
    """ Base observer class that is called in the simulation by the particle class to
    determine when to write out (observe) data. 
     
    The conditions to observe can be based on the time (or step) or the coordinates of 
    the particle.
    - step number [unique_steps] -> time (TimeEvolutionObservers)
    - radius of observer sphere [shperes] -> sphere around source (SphericalObservers)
    - cartesian coordinates [box_dimensions] -> box around source (BoxObserver) (not yet implemented)
    All special observer will create an Observer object and specify the relevant parameters
    for the observation conditions (unique_steps, shperes, box_dimensions)

    Attributes:
        substeps: An b1 array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        all_steps: A bool that determines if all steps should be observed, independet of steps array.
        steps: A float32 array specifying all steps that should be observed.
        pos: A float32 array for the position of the current particle.
        spheres: A float32 array for specifying the radii of the observer spheres.
        ps: ParticleState of the current particle.
 
    """

    def __init__(self, steps, substeps, spheres):
        self.substeps = substeps
        self.spheres = spheres
        self.all_steps = False
        if -1 in steps:
            # the -1 is the key to say that all steps should be observed
            self.all_steps = True
        self.steps = self.get_unique_steps(steps)
        print('Observer initialized')
    

    def get_unique_steps(self, steps):
        """Remove all step number duplicates. 

        Special observer such as the TimeEvolutionObserverLog() may introduce duplicate 
        step numbers [1,1,1,...] that lead to a wrong number of steps during the description 
        output because the step 1 will be only observed once, while it is counted several 
        times in the description output.

        Args:
            steps (int32[:]): All step numbers that should be observed.
        Returns:
            unique_steps (int32[:]): steps but removed all duplicates.
        """

        unique_steps = []
        for s in steps:
            if s not in unique_steps:
                unique_steps.append(s)
        return np.array(unique_steps, dtype=np.int32)
        

    def observe(self, ps):
        """Check if and what to observe.
        
        Decides if the current particle state should be observed based on the criterions specified 
        in the observer instance.

        Args:
            ps (ParticleState): Current particle state.
        Returns:
            new data row or None
        """

        if ps.substep == 2 and len(self.spheres) > 1:
            data_on_sphere = self.check_on_sphere(ps)
            if data_on_sphere != None:
                return data_on_sphere
        if self.substeps[ps.substep]:
            if self.all_steps or ps.step in self.steps:
                return self.data_row(ps, -1.0)
            else:
                return None
        else:
            return None


    def check_on_sphere(self, ps):
        """Check if particle passes a sphere during the propagation step.
        
        Checks if the particle crosses through a sphere with radii r in the current
        propagation step. As there may be many spheres, a loop over all user-
        specified radii is needed.

        Args:
            ps (ParticleState): Current particle state.
        Returns:
            new data row or None
        """
        radius_2 = 0
        radius_prev_2 = 0
        for i in range(ps.dimensions):
            radius_2 = radius_2 + ps.pos[i]**2
            radius_prev_2 = radius_prev_2 + ps.pos_prev[i]**2
        radius = radius_2**0.5
        radius_prev = radius_prev_2**0.5
        for j in range(1, len(self.spheres)):
            r_s = self.spheres[j]
            if (radius > r_s and radius_prev < r_s) or (radius < r_s and radius_prev > r_s):
                return self.data_row(ps, r_s)
        
        return None


    def data_row(self, ps, radius):
        """List of particle parameters that should be observed.

        Args:
            ps (ParticleState): Current particle state.
            radius (float32): Radius of the sphere that the particle is located at.
        Returns:
            data_row_list (list): All items that should be observed. 
        """

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
        """Description including that of the parameters.
        
        Note: description does not contain the information of the underling special observer 
        (if there was one) that was used during the observer initialization. To get 
        this info, get_description(self) has to be called directly on the instance of 
        the special observer class (see tutorials for details).
        """
        self.get_description_general()
        self.get_description_parameters()


    def get_description_general(self):
        """General description.

        Called by all special observer classes below.
        Introduction of the description output.
        """
        print("""Description Observer:
                The observer defines the conditions for when to write data to the output.\n""")


    def get_description_parameters(self):   
        """Paremeters description including the values set.

        Called by all special observer classes below.
        Introduction of the description output.
        """
        print('spheres: ' , self.spheres)
        print('steps [0:10]: ' , self.steps[0:10])
        print('steps [-11:-1]: ' ,self.steps[-11:-1])
        print('nr steps: ' , len(self.steps))
        print('substeps: ', self.substeps)  
        print('all_steps: ', self.all_steps)         


#-----------------------------------------------------------------------------
"""Special observer classes

Below are the abstract base class and all sub classes of the special observers
that have to be added to the simulation. Each special observer stores a
Observer object in its instance parameter to be used in the simulation.
This diversions is needed because numba does not support 
inheritance via ABC and Propagator() needs the label @jitclass as it is called 
during the numba optimized simulation loop of the run_simulation() function. 
This workaround supports both concepts with the 
advantages of fast code and easy addition of new observers where the structure 
is now defined by the Abstract Base class and enforeced via the ABCMeta class
"""


class AbstractSpecialObserverMeta(ABCMeta):
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
        obj = super(AbstractSpecialObserverMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj



class AbstractSpecialObserver(object, metaclass=AbstractSpecialObserverMeta):
    """Abstract base class for all special observers.
    
    Functions with the label @abstractmethod have to be implemented in the special 
    observer classes.

    Attributes:
        substeps_bool: An b1 array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        steps: A float32 array specifying all steps that should be observed.
        spheres: A float32 array for specifying the radii of the observer spheres.
        column: A float32 array for specifying the column names of the data output.
        observer: The special observer.
    """
    
    required_attributes = [
        'steps', 
        'substeps_bool',
        'spheres'
    ]
 
    @abstractmethod
    def __init__(self, order):
        """ Implementation required in all sub classes. All required_attributes have to 
        be implemented in sub classes.
        """
        pass

    def init_observer(self, substeps, spheres):
        """Initialize the Observer 
        
        Sets the important parameters and calls the @abstractmethods that are implemented
        in each special observer class that are derived from the current abstract base class.

        Args:
            substeps (float32[:]): Specifying which substeps need to be observed.
            spheres (float32[:]): Specifying the radii of the observer spheres.
        """

        self.column = ['id', 'i', 'd', 'x', 'y', 'z', 'phi', 'pitch_angle', 'radius', 'sub_step']
        self.substeps_bool = np.array(substeps)
        self.steps = self.set_steps_int() 
        self.spheres = np.array(spheres, dtype=np.float32)
        # have to store all relevant observation parameters in the Observer class that 
        # has the @jitclass label from numba. This is important, as the Particle class is also 
        # labeled with @jitclass and can thus only call @jitclass classes. The usage of numba is 
        # for performance resaons.
        self.observer = Observer(self.steps, self.substeps_bool, self.spheres)

    @abstractmethod
    def set_steps(self):
        """Sets the number of steps that should be observed.
        """
        pass

    def set_steps_int(self):
        """Converts list of steps to np array of ints that is known in the @jitclass Observer().
        """
        steps = self.set_steps()
        steps_int32 = np.array(steps, dtype=np.int32)
        return steps_int32
    
    @abstractmethod
    def get_description_observer_type(self):
        """Gives the name of each special observer.
        """
        pass

    def get_description(self):
        """Prints a description of the observer with the status of all attributes.
        """
        self.observer.get_description_general()
        self.get_description_observer_type()
        self.observer.get_description_parameters()

    def get_column_names(self):
        """Gets called in the simulation script and returns the column names of the data output.
        """
        return self.column

    

class ObserverAllSteps(AbstractSpecialObserver):
    """Observes particles in all propagation steps.

    Attributes:
        substeps: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        steps_input: A float array specifying all steps that should be observed.
    """

    def __init__(self, substeps):
        spheres = [-1.0]
        self.init_observer(substeps, spheres)

    def set_steps(self):
        steps = [-1]
        return steps

    def get_description_observer_type(self):
        print('observer tpye: ObserverAllSteps')



class TimeEvolutionObserverLog(AbstractSpecialObserver):
    """Observes particles at the user specified step numbers.

    The user only gives the minimum, the maximum and the total step numbers. The
    TimeEvolutionObserverLog computes the list (logarithmically).

    Attributes:
        substeps: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        min_steps: An int that gives the minimal step number that should be observed.
        max_steps: An int that gives the maximal step number that should be observed.
        nr_steps: An int that gives the number of steps that should be observed.
        steps_input: A float array specifying all steps that should be observed.
    """

    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        self.min_steps = min_steps
        self.max_steps = max_steps
        self.nr_steps = nr_steps
        spheres = [-1.0]

        self.init_observer(substeps, spheres)
        
    def set_steps(self):
        """Sets the number of steps that should be observed.

        Takes the minimum and maximum step numbers into account to generate the list of
        step numbers based on the given number of steps. Here, the steps are spaced 
        logarithmically.
        """
        steps = np.logspace(np.log10(self.min_steps), np.log10(self.max_steps), self.nr_steps)
        return steps

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserverLog')



class TimeEvolutionObserverLin(AbstractSpecialObserver):
    """Observes particles at the user specified step numbers.

    The user only gives the minimum, the maximum and the total step numbers. The
    TimeEvolutionObserverLin computes the list (linearly).

    Attributes:
        substeps: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        min_steps: An int that gives the minimal step number that should be observed.
        max_steps: An int that gives the maximal step number that should be observed.
        nr_steps: An int that gives the number of steps that should be observed.
        steps_input: A float array specifying all steps that should be observed.
    """

    def __init__(self, min_steps, max_steps, nr_steps, substeps):
        self.min_steps = min_steps
        self.max_steps = max_steps
        self.nr_steps = nr_steps
        spheres = [-1.0]

        self.init_observer(substeps, spheres)

    def set_steps(self):
        """Sets the number of steps that should be observed.

        Takes the minimum and maximum step numbers into account to generate the list of
        step numbers based on the given number of steps. Here, the steps are spaced 
        linearly.
        """
        steps = np.linspace(self.min_steps, self.max_steps, self.nr_steps)
        return steps

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserverLin')
    


class TimeEvolutionObserver(AbstractSpecialObserver):
    """Observes particles at the user specified step numbers.

    The user passes the list of steps to the TimeEvolutionObserver.

    Attributes:
        substeps: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        steps_input: A float array specifying all steps that should be observed.
    """

    def __init__(self, steps_input, substeps):
        self.steps_input = steps_input
        spheres = [-1.0]

        self.init_observer(substeps, spheres)

    def set_steps(self):
        """Sets the number of steps that should be observed.
        """
        steps = self.steps_input 
        return steps

    def get_description_observer_type(self):
        print('observer tpye: TimeEvolutionObserver')



class SphericalObserver(AbstractSpecialObserver):
    """Observes particles on a spheres with given radii.
    
    When particles pass through the sphere, they will be observed.

    Attributes:
        substeps: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        spheres: A float array for specifying the radii of the observer spheres.
    """

    def __init__(self, substeps, spheres):
        self.steps_input = []
        self.spheres = [-1.0] + spheres

        self.init_observer(substeps, self.spheres)

    def get_description_observer_type(self):
        print('observer tpye: SphericalObserver')