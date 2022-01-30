"""The Source initializes the particles in the beginning of the simulation.

Background information:
    Current simulation codes propagate charged particles through magnetic fields 
    either by ballistic codes solving the Lorentz force at each step or by using 
    particle distributions arising from the solution of the classical diffusion 
    equation. Whereas the former method provides an accurate description of the 
    particle trajectories using small step sizes, this method is computationally 
    expensive, especially for low particle energies. The method based on the 
    diffusion equation, on the other hand, is extremely fast since the particle 
    distribution is analytically given at each time point, but by the nature of 
    the method, particles can only be described statistically. Even using quasi-
    particles, the individual trajectories are useless. For applications in which 
    statistical statements and averaging over many particles are intended, this 
    method is preferred in many areas of astronomy because of the short 
    simulation times.
    It is important to note, however, that the diffusion equation guarantees an 
    adequate description of the particles only in the limit of infinitely large 
    times. Numerically, however, this approach can also be used starting from 
    times for which a diffusive behavior occurs. This typically happens after a 
    simulation time which is in the order of magnitude of the mean free path 
    length. 
    Consequently, the use of propagation codes based on the diffusion equation 
    is not suitable in compact objects. Examples include AGNs, supernovas, 
    pulsars, etc.

Explanation:
    In analogy to existing diffusion codes, the following routine requires the 
    information of the diffusion coefficients. Here, however, we do not use the 
    solution of the diffusion equation, since this only provides a particle 
    distribution and makes statistical statements about individual particle 
    trajectories. Instead, we propagate particles according to the two-step 
    propagation routine. For strong turbulence levels we 
    use the correlated random walk (CRW) in Cartesian coordinates and for weak 
    turbulence levels the CRW in cylindrical coordinates.

Typical usage example:

    import proppy as pp
    
    nr_steps = 1*10**4
    step_size = 0.5*10**10 # [m]
    mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
    rms = 1 # Gaus
    magnetic_field = pp.OrderedBackgroundField(rms, [0,0,1]).magnetic_field

    propagator = pp.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    
    sim = pp.Simulation()
    sim.add_propagator(propagator)
    sim.propagator.get_description()
"""


from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from .magnetic_field import *
from .particle_state import *
from .constants import *
from abc import ABCMeta, abstractmethod


propagation_spec = [
    ('cartesian', b1),
    ('cylindrical', b1),
    ('nr_steps', int32),
    ('dimensions', int32),
    ('background_direction', int32),
    ('step_distance', float32),
    ('step_size', float32),
    ('step_size_diff_factor', float32),
    ('gyroradius_eff', float32),
    ('speed', float32),
    ('pitch_angle_const', b1),

    ('prob', float32[:]),

    ('pos', float32[:]),
    ('direction', float32[:]),
    ('mfp', float32[:]),
    ('phi', float32),
    ('pitch_angle', float32),
    ('distance', float32),

    ('magnetic_field', MagneticField.class_type.instance_type),
    ('constants', Constants.class_type.instance_type),
    ('particle_state', ParticleState.class_type.instance_type),
]

@jitclass(propagation_spec)
class Propagator():
    """ Base propagator class that is called in the simulation by the particle class to
    determine the next particle state. 
     
    In analogy to existing diffusion codes, the following routine requires the 
    information of the diffusion coefficients. Here, however, we do not use the 
    solution of the diffusion equation, since this only provides a particle 
    distribution and makes statistical statements about individual particle 
    trajectories. Instead, we propagate particles according to the two-step 
    propagation routine. 
    - For strong turbulence levels: We use the correlated random walk (CRW) 
    in Cartesian coordinates.
    - For weak turbulence levels: We use the CRW in cylindrical coordinates.

    Attributes: (for types, see the propagation_spec above)
        speed: Speed of particle in [m/s].
        cartesian: Should use Cartesian coordinates?
        cartesian: Should use cylindrical coordinates?
        nr_steps: Number of steps for each particle in simulation.
        step_size: Size of each individual step in [m].
        dimensions: Number of dimensions.
        pitch_angle_const: Should the pitch angle remain const or should there be pitch angle scattering?
        background_direction: Direction of a background field (0=x-axis,...).
        magnetic_field: The magnetic field.
        prob: Probability to change in each direction in a step.
        mfp: Mean free paths of particles in each direction in [m].
        step_size_diff_factor: Increase step size in diffusive phase by this factor. 
    """

    def __init__(self, nr_steps, step_size, magnetic_field, cartesian, mfp, step_size_diff_factor=1.0, constants=Constants()):
        print('Propagator initialized')
        self.cartesian = cartesian
        self.cylindrical = not cartesian
        self.nr_steps = nr_steps
        self.step_size = step_size # [m]
        self.dimensions = 3
        self.pitch_angle_const = True
        self.background_direction = 2
        self.magnetic_field = magnetic_field
        self.mfp = mfp
        self.step_size_diff_factor = step_size_diff_factor
        self.constants = constants
        self.prob = self.set_prob_init(self.mfp, self.constants.speed, self.step_size)

    
    def set_prob_init(self, mean_free_path, speed, step_size):
        """Calculate the propabilities to change directions.

        The propabilities to change directions are based on the mean free paths
        that depend on the diffusion coefficients.

        Args:
            mean_free_path: Mean free paths of particles in each direction in [m].
            speed: Speed of the particles in [m/s].
            step_size: Size of the steps in [m].

        Returns:
            probabilities: Probability to change the direction in one prop. step.
        """
        xi = [
            speed / mean_free_path[0] / 2.0, 
            speed / mean_free_path[1] / 2.0, 
            speed / mean_free_path[2] / 2.0
        ] # [1/s] frequency of change
        tau_step = step_size / speed
        propabilities = [xi[0] * tau_step, xi[1] * tau_step, xi[2] * tau_step]
        return np.array(propabilities, dtype=np.float32)
        

    def change_direction(self, direction, particle_state):
        """Change particle direction during step in random walk.
        
        Change in direction happens with a propability that is defined by the 
        diffusion coeffcients (used to compute the mean free paths). The 
        information about the probabilities is stored in self.prob.

        Args:
            direction: Current direction of particle.

        Returns:
            direction: New direction of particle that is the same as the input direction (-> correlated random walk).
        """
        for p in range(self.dimensions):
            prop = self.get_adaptive_prob(particle_state, self.prob[p], p)
            if np.random.random() < prop:
                direction[p] = -1*direction[p]
        return direction


    def change_pitch_angle(self, pitch_angle):
        """Change pitch angle of particles.

        Only used when pitch angle shouldn't kept constant! Changes in the pitch angle are 
        caused by resonant scattering of particles at fluctuations of the turbulence that 
        satisfy the resonance scattering criterion. These changes in pitch angle are 
        approximated in Kulsrud & Pearce (1969, ApJ, 156, 445) and 
        Reichherzer et al. (2020, MNRAS) as follows:
        delta mu = b/B.
        Here, b is the rms field strength of the turbulence and B the magnetic field strength 
        of the ordered magnetic field lines. Only valid for weak turbulence levels b << B.
        
        Args:
            pitch_angle: Current pitch angle of particle.

        Returns:
            new_pitch_angle: New pitch angle caused by pitch angle scattering.
        """
        if self.pitch_angle_const or np.random.random() >= self.prob[self.background_direction]:
            # if the pitch angle should be constant no calculations needed here
            return pitch_angle
        # mu = cos(pitch angle) -> delta mu = cos(theta_1) - cos(theta_0)
        delta_mu = 0.1 # 0.1 corressponds to a weak turbulence level b/B = 0.1
        if np.random.random() < 0.5:
            # change sign of delta_mu with probability of 50%. Adding and subtraction from the pitch_angle 
            # has the same probability
            delta_mu = -delta_mu
        # delta mu = cos(theta_1) - cos(theta_0)
        # -> theta_1 = arccos(delta_mu - cos(theta_1))
        pitch_angle_0 = pitch_angle
        pitch_angle_1 = np.arccos(delta_mu + np.cos(pitch_angle_0))
        new_pitch_angle = pitch_angle_1
        return new_pitch_angle


    def move_substep(self, particle_state):
        """Move one substep.
        
        The number of the current substep is stored in the particle state.
        
        Args:
            particle_state: Current particle state.
            
        Returns: 
            new_particles_state: New particle state after propagation of substep.
        """
        # local step
        particle_state, move_local_array = self.move_local(particle_state)
        # global step
        new_particle_state = self.move_global(particle_state, move_local_array)
        return new_particle_state


    def move_local(self, particle_state):
        """Movements during a substep in the local frame.
        
        The local frame is defined such that the magnetic field points into the
        z-axis. Depending on the variables self.cartesian and self.cylindrical, 
        the particle moves via a correlated random walk in Carteesian or 
        cylindrical coordinates.
        
        Args:
            particle_state: Current particle state.
            
        Returns: 
            new_particles_state: New particle state after propagation of substep in local frame.
        """
        if self.cartesian:
            # cartesian coordinates -> move in x, y and z directions
            particle_state, move_local = self.move_cartesian_simple(particle_state)
        else:
            # cylindrical coordinates -> move in phi, rho and z directions
            if particle_state.substep == 0:
                particle_state, move_local = self.move_phi(particle_state)
            if particle_state.substep == 1:
                particle_state, move_local = self.move_rho(particle_state)
            if particle_state.substep == 2:
                particle_state, move_local = self.move_cartesian(particle_state)
        return particle_state, move_local


    def move_global(self, ps, move_local):
        """Movements during a substep in the global frame -> new position.
        
        Currently only magnetic fields in z-axis are supported. Therfore, the
        movement in the local frame is the same as the movement in the global
        frame. Therefore, the new position can be derived by the local movement
        without any transformations.
        
        Args:
            ps: Current particle state.
            move_local: An array that describes the local move in the substep.
            
        Returns: 
            ps: New particle state after propagation of substep in global frame.
        """
        for s in range(self.dimensions):
            ps.pos[s] = ps.pos[s] + move_local[s]
        return ps
    
                 
    def move_cartesian_simple(self, particle_state):
        """Movement in substep in Cartesian coords.
        
        During the complete propagation step, the particle moves one step size further, without taking the pitch angle into account
        
        Args:
            particle_state: Current particle state.
            
        Returns: 
            particles_state: New particle state after propagation of substep in global frame.
            move_local: An array that describes the local move in the substep.
        """
        step_size = self.get_adaptive_step_size(particle_state)
        distance_s = particle_state.direction[particle_state.substep] * step_size / 3**0.5
        particle_state.distance += step_size / 3.0
        move_local = [0,0,0]
        s = particle_state.substep
        move_local[s] = distance_s
        return particle_state, self.float_array(move_local)


    def get_adaptive_step_size(self, particle_state):
        if particle_state.distance > 5*self.mfp[particle_state.substep]:
            adaptive_step_size = self.step_size * self.step_size_diff_factor
        else:
            adaptive_step_size = self.step_size
        return adaptive_step_size


    def get_adaptive_prob(self, particle_state, prob, p):
        if particle_state.distance > 5*self.mfp[p]:
            adaptive_prob = prob * self.step_size_diff_factor
        else:
            adaptive_prob = prob
        return adaptive_prob


    
    def move_cartesian(self, particle_state):
        """Movement in substep in Cartesian coords. or in z-axis in cylindrical coords.
        
        During the complete propagation step, the particle moves one step size further.
        In components in each axis depends on the pitch angle of the particle.
        
        Args:
            particle_state: Current particle state.
            
        Returns: 
            particles_state: New particle state after propagation of substep in global frame.
            move_local: An array that describes the local move in the substep.
        """
        distance_s = 0.0
        if particle_state.substep == self.background_direction:
            distance_s = self.step_size * np.cos(particle_state.pitch_angle) * particle_state.direction[particle_state.substep]
        else:
            distance_s = self.step_size * np.sin(particle_state.pitch_angle) / 2**0.5 * particle_state.direction[particle_state.substep]
        if self.cartesian:
            particle_state.distance = particle_state.distance + np.abs(distance_s)
        else:
            particle_state.distance = particle_state.distance + self.step_size
        move_local = [0,0,0]
        s = particle_state.substep
        move_local[s] = distance_s
        return particle_state, self.float_array(move_local)
        
        
    def move_phi(self, particle_state):
        """Movement in substep along the phi-direction in cylindrical coords.
        
        Here, the particles phi angle is increased or decreased as defined in the current
        direction array for the phi angle. The particle moves in positive or negative 
        direction in the circle with a radius that is defined by the effective gyroradius.
        
        Args:
            particle_state: Current particle state.
            
        Returns: 
            particles_state: New particle state after propagation of substep in global frame.
            move_local: An array that describes the local move in the substep.
        """
        phi_old = particle_state.phi
        delta_phi = self.compute_delta_phi(particle_state)
        particle_state.phi = phi_old + delta_phi * particle_state.direction[0]
        chi_x_1 = particle_state.gyroradius_eff * (np.cos(particle_state.phi) - np.cos(phi_old))
        chi_y_1 = particle_state.gyroradius_eff * (np.sin(particle_state.phi) - np.sin(phi_old))
        move_local = [chi_x_1, chi_y_1, 0]
        return particle_state, self.float_array(move_local)

                      
    def move_rho(self, particle_state):
        """Movement in substep along the rho-direction in cylindrical coords.
        
        Here, the particle move in the rho direction that is specified by the current
        direction array.
        
        Args:
            particle_state: Current particle state.
            
        Returns: 
            particles_state: New particle state after propagation of substep in global frame.
            move_local: An array that describes the local move in the substep.
        """
        #particle_state.distance = particle_state.distance + np.abs(distance_s)
        delta_rho = self.step_size * np.sin(particle_state.pitch_angle) / 2**0.5
        chi_x_2 = np.cos(particle_state.phi) * particle_state.direction[1] * delta_rho
        chi_y_2 = np.sin(particle_state.phi) * particle_state.direction[1] * delta_rho
        move_local = [chi_x_2, chi_y_2, 0]
        return particle_state, self.float_array(move_local)


    def compute_delta_phi(self, ps):
        """Compute the change in the phi angle that is needed based on the pitch angle.
        
        Args:
            particle_state: Current particle state.
            
        Returns: 
            delta_phi: Change in phi angle needed during the phi substep in cylindrical coords.
        """
        delta_rho = self.step_size * np.sin(ps.pitch_angle)
        delta_phi = 2 * np.arcsin(delta_rho / (2 * 2**0.5 * ps.gyroradius_eff))
        return delta_phi


    def float_array(self, pos):
        return np.array([pos[0], pos[1], pos[2]], dtype=np.float32)


    def set_constants(self, constants):
        self.constants = constants


    def set_gyroradius(self, ps):
        """Compute gyroradius given the magnetic field strength and the particle energy.
        
        It is important to notice that we have to consider an effective gyroradius for the
        movement of the particle in phi direction, because combining the movement in phi and
        rho directions gives an effective gyromation with a larger gyroradius than the one
        used or the phi movement. In order to have a gyromotion with the correct gyroradius
        of the combined substeps in rho and phi direction, the gyroradius used for the phi
        movement has to be smaller -> effective gyrorradius.
        Correcting for moving in rho direction (perp to phi) --> gyration increases by 2**0.5,
        which is why we have to divide here. 

        Args:
            ps: Current particle state.
            
        Returns: 
            ps: New particle state with an effective gyroradius for the current energy and b-field strength.
        """
        # Normalization factor for gyroradius for protons (v=c) with 1eV in magnetic 
        # field with strength 1Gaus.
        gyroradius_0 = 3.336*10**(-5) # in [m]
        if self.magnetic_field.rms == 0:
            ps.gyroradius = gyroradius_0
        else:
            ps.gyroradius = gyroradius_0 * ps.energy / self.magnetic_field.rms # in [m]
        ps.gyroradius_eff = ps.gyroradius / 3**0.5 
        return ps
        

    def set_pitch_angle_const(self, const_bool):
        """Keep the pitch angle either constant or allow for changes 
        during each propagation step.
        """
        self.pitch_angle_const = const_bool


    def set_dimensions(self, dimensions):
        """Default is 3d -> dimensions = 3. More than 3 dimensions are not supported.
        """
        self.dimensions = dimensions


    def set_cartesian_coords(self, cartesian_bool):
        """Choose Cartesian coords or not.
        
        There are cartesian or cylindrical coordinates available. 
        cylindrical coordinates activated lead to the usage of
        move_phi, move_rho and move_cartesian for the z-direction
        """
        self.cartesian = cartesian_bool
        self.cylindrical = not cartesian_bool

    
    def set_cylindrical_coords(self, cylindrical):
        """Choose cylindrical coords or not.

        There are cartesian or cylindrical coordinates available. 
        Cylindrical coordinates activated lead to the usage of
        move_phi, move_rho and move_cartesian for the z-direction
        """
        self.cartesian = not cylindrical
        self.cylindrical = cylindrical


    def set_nr_steps(self, nr_steps):
        """Change number of steps.
        """
        self.nr_steps = nr_steps

    
    def set_step_size(self, step_size):
        """Change distance of each step that particles travel.
        
        units = [m]
        """
        self.nr_steps = step_size


    def set_prob(self, mean_free_path):
        self.prob = self.set_prob_init(mean_free_path, self.constants.speed, self.propagator.step_size)
    

    def set_magnetic_field(self, magnetic_field):
        self.magnetic_field = magnetic_field


    def get_description(self):
        """Description.
        
        Note: description does not contain the information of the underling special propagator 
        (if there was one) that was used during adding the propagator to the simulation. To get 
        this info, get_description(self) has to be called directly on the instance of the used 
        propagator (see tutorials for details).
        """
        self.get_description_general()
        self.get_description_parameters()


    def get_description_general(self):
        """Print out the discription of the object with all relevant instance parameters.
        """
        print("""Description Propagator:
                The propagator is responsible for the movement of the particles. 
                It performs the change of direction and the movement in the respective direction.
                There are two phases:
                 - change direction with probability (see below)
                 - move in all directions
                The movement takes place according to the correlated random walk (CRW).\n""")


    def get_description_parameters(self):
        """Print out the discription of the object with all relevant instance parameters.
        """
        print("""Description Propagator:
                The propagator is responsible for the movement of the particles. 
                It performs the change of direction and the movement in the respective direction.
                There are two phases:
                 - change direction with probability (see below)
                 - move in all directions
                The movement takes place according to the correlated random walk (CRW).\n""")
        if self.cartesian:
            print('coordinates: Cartesian coordinates')
        else:
            print('coordinates: Cylindrical coordinates')
        print('dimensions: ', self.dimensions)
        if self.dimensions < 3:
            print(self.dimensions, ' dimensions not supoorted yet!')
        if self.pitch_angle_const:
            print('pitch angle: constant')
        else:
            print('pitch angle: not constant')
        print(self.constants.print_unit_speed())
        print('speed: ', self.constants.speed, self.constants.print_unit_speed())
        print('number steps: ', self.nr_steps)  
        print('step size: ', self.step_size, self.constants.print_unit_distance())  
        print('step duration: ', self.step_size / self.constants.speed, ' s') 
        print('total distance: ', self.step_size * self.nr_steps, self.constants.print_unit_distance())
        print('total duration: ', self.step_size * self.nr_steps / self.constants.speed, ' s')
        print('probability to change directions in step: ', self.prob*100, '%')  




#---------------------------------------------------------------------------------
"""Below are the abstract base class and all sub classes of the special propagators
that have to be added to the simulation. Each special propagator stores a
Propagator object in its instance parameter to be used in the simulation.This 
diversions is needed because numba does not support inheritance via ABC and 
Propagator() needs the label @jitclass as it is called during the numba optimized 
simulation loop of the run_simulation() function. This workaround supports both 
concepts with the advantages of fast code and easy addition of new propagator where 
the structure is now defined by the Abstract Base class and enforeced via the 
ABCMeta class.
"""



class AbstractPropagatorMeta(ABCMeta):
    """Abstract meta class to check if all required attributes are implemented in the 
    sub classes.
    """
    required_attributes = []

    def __call__(self, *args, **kwargs):
        """ Checks if required attributes that have to be implemented in __init__ of all
        sub classes are really implemented. 

        Raises:
            ValueError: an error if not all required attributes are implemented.
        """
        obj = super(AbstractPropagatorMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj



class AbstractPropagator(object, metaclass=AbstractPropagatorMeta):
    """Abstract base class for all special propagator.
    
    Functions with the label @abstractmethod have to be implemented in the special 
    propagator classes.

    Attributes: 
        propagator: The special propagator.
        dimensions: Number of dimensions.
        speed: Speed of the particles in [m/s].
        mfp: Mean free paths of particles in each direction in [m].
        nr_steps: Number of steps for each particle in simulation.
        magnetic_field: The magnetic field.
        step_size: Size of the steps in [m].
        isotropic_diffusion: Is the particle diffusion isotropic?
    """
    
    required_attributes = [
        'propagator', 
        'dimensions', 
        'speed', 
        'mfp', 
        'nr_steps', 
        'magnetic_field', 
        'step_size', 
        'isotropic_diffusion',
        'cartesian',
        'cylindrical'
    ]
 
    @abstractmethod
    def __init__(self, order):
        """Implementation required in all sub classes. All required_attributes have to 
        be implemented in sub classes.
        """
        pass

    def set_basic_parameters(self):
        self.dimensions = 3
        self.speed = 2.998*10**8 # [m/s]


    def set_pitch_angle_const(self, const_bool):
        """Keep the pitch angle either constant or allow for changes during each propagation step.
        """
        self.pitch_angle_const = const_bool
        self.propagator.pitch_angle_const = const_bool


    def set_pitch_angle_const(self, const_bool):
        """Keep the pitch angle either constant or allow for changes 
        during each propagation step.
        """
        self.pitch_angle_const = const_bool


    def set_dimensions(self, dimensions):
        """Default is 3d -> dimensions = 3. More than 3 dimensions are not supported.
        """
        self.dimensions = dimensions
        self.propagator.dimensions = dimensions        


    def set_cartesian_coords(self, cartesian_bool):
        """Choose Cartesian coords or not.
        
        There are cartesian or cylindrical coordinates available. 
        cylindrical coordinates activated lead to the usage of
        move_phi, move_rho and move_cartesian for the z-direction
        """
        self.cartesian = cartesian_bool
        self.cylindrical = not cartesian_bool
        self.propagator.cartesian = cartesian_bool
        self.propagator.cylindrical = not cartesian_bool

    
    def set_cylindrical_coords(self, cylindrical):
        """Choose cylindrical coords or not.

        There are cartesian or cylindrical coordinates available. 
        Cylindrical coordinates activated lead to the usage of
        move_phi, move_rho and move_cartesian for the z-direction
        """
        self.cartesian = not cylindrical
        self.cylindrical = cylindrical
        self.propagator.cartesian = not cylindrical
        self.propagator.cylindrical = cylindrical


    def set_speed(self, speed):
        """Set particle speed.
        
        units = [m/s]
        change the speed of the particles. The default speed is the speed of light 
        that is valid for relativistic particles.
        """
        self.speed = speed
        self.propagator.speed = speed


    def set_nr_steps(self, nr_steps):
        """Change number of steps.
        """
        self.nr_steps = nr_steps
        self.propagator.nr_steps = nr_steps

    
    def set_step_size(self, step_size):
        """Change distance of each step that particles travel.
        
        units = [m]
        """
        self.nr_steps = step_size
        self.propagator.nr_steps = step_size        


    def set_prob(self, mean_free_path):
        self.propagator.prob = self.set_prob_init(mean_free_path, self.propagator.constants.speed, self.propagator.step_size)
        self.prob = self.propagator.prob
    

    def set_magnetic_field(self, magnetic_field):
        """Allow to set the propagator magnetic field.
        """
        self.magnetic_field = magnetic_field
        self.propagator.magnetic_field = magnetic_field


    def convert_mfp_input(self, mfp_input):
        """Check the input of the mean free paths.
        """
        if self.isotropic_diffusion == False:
            if isinstance(mfp_input, float) or isinstance(mfp_input, int) or len(mfp_input) < self.dimensions:
                # error handeling by wrong input
                raise ValueError('Input error: please provide a list of anisotropic mean free paths as a list with length of dimensions (default length: 3).')
            else:
                return mfp_input
        mfp = []
        if isinstance(mfp_input, float) or isinstance(mfp_input, int):
            for i in range(self.dimensions):
                mfp.append(mfp_input)
        else:
            for i in range(self.dimensions):
                mfp.append(mfp_input[0])

        return np.array(mfp, dtype=np.float32)

    
    def init_jitclass_propagator(self):
        """Initialize the Propagator. 
        
        Sets the important parameters and calls the @abstractmethods that are implemented
        in each special propagator class that are derived from the current abstract base class.
        """
        self.set_basic_parameters()
        self.mfp = self.convert_mfp_input(self.mfp)
        propagator = Propagator(self.nr_steps, self.step_size, self.magnetic_field, self.cartesian, self.mfp)
        # have to store all relevant propagation parameters in the Propagator class that 
        # has the @jitclass label from numba. This is important, as the Particle class is also 
        # labeled with @jitclass and can thus only call @jitclass classes. The usage of numba is 
        # for performance resaons.
        self.propagator = propagator

    
    @abstractmethod
    def get_description_propagator_type(sefl):
        pass


    def get_description(self):
        """Description.
        
        Print the information of the relevant parameters and the description of 
        the special propagation type that was chosen.
        """
        self.propagator.get_description_general()
        self.get_description_propagator_type()
        self.propagator.get_description_parameters()


    def get_description_general(self):
        """Called by all special propagator classes below. Introduction of the description output.
        """
        print("""Description Propagator:
                The propagator is responsible for the movement of the particles. 
                It performs the change of direction and the movement in the respective direction.
                There are two phases:
                 - change direction with probability (see below)
                 - move in all directions
                The movement takes place according to the correlated random walk (CRW).\n""")


    def get_description_parameters(self):   
        """Called by all special propagator classes below. Print out all relevant 
        instance parameters.
        """
        print('number steps: ', self.nr_steps)
        print('call get_description directly on the propagator that was added to the simulation:\n')
        print('sim = proppy.Simulation()')
        print('...')
        print('sim.add_propagator(propagator)')
        print('sim.propagator.get_description()')



class IsotropicPropagatorDefault(AbstractPropagator):
    """Default isotropic propagator.

    Usage of no background magnetic field with only a turbulent one. Isotropic
    diffusion coefficients with corresponding mfp = 10**12 m. Propagating with
    Cartesian coordinates.

    Attributes:
        nr_steps: Number of steps for each particle in simulation.
        step_size: Size of the steps in [m].
        mfp: Mean free paths of particles in each direction in [m].
        magnetic_field: The magnetic field.
        isotropic_diffusion: Is the particle diffusion isotropic?
    """

    def __init__(self):
        self.nr_steps = 2*10**5
        self.step_size = 0.5*10**10 # [m]
        # isotropic diffusion coefficient
        self.mfp = np.array([10**12, 10**12, 10**12], dtype=np.float32)  # [m]
        # no background magnetic field
        rms = 0
        self.magnetic_field = DefaultBackgroundField(rms).magnetic_field
        self.isotropic_diffusion = True
        self.cartesian = True
        self.cylindrical = False

        self.init_jitclass_propagator() 


    def get_description_propagator_type(self):
        print('propagation tpye: IsotropicPropagatorDefault')



class IsotropicPropagator(AbstractPropagator):
    """Isotropic propagator with user-specified parameters.

    Usage of no background magnetic field with only a turbulent one. Isotropic
    diffusion coefficients and therefore isotropic mfp. Propagating with
    Cartesian coordinates. User can specify mfp, nr steps and step size.

    Attributes:
        nr_steps: Number of steps for each particle in simulation.
        step_size: Size of the steps in [m].
        mfp: Mean free paths of particles in each direction in [m].
        magnetic_field: The magnetic field.
        isotropic_diffusion: Is the particle diffusion isotropic?
    """

    def __init__(self, mfp, nr_steps, step_size, step_size_diff_factor=1.0):
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.mfp = mfp
        self.step_size_diff_factor = step_size_diff_factor
        # no background magnetic field
        rms = 0
        self.magnetic_field = DefaultBackgroundField(rms).magnetic_field
        self.isotropic_diffusion = True
        self.cartesian = True
        self.cylindrical = False

        self.init_jitclass_propagator()


    def get_description_propagator_type(self):
        print('propagation tpye: IsotropicPropagator')
  


class AnisotropicPropagator(AbstractPropagator):
    """Anisotropic propagator with user-specified parameters.

    Usage of a background magnetic field and a turbulent one. Anisotropic
    diffusion coefficients and therefore anisotropic mfp. Propagating with
    Cylindircal coordinates. User can specify mfp, nr steps and step size,
    as well as magnetic field.

    Attributes:
        nr_steps: Number of steps for each particle in simulation.
        step_size: Size of the steps in [m].
        mfp: Mean free paths of particles in each direction in [m].
        magnetic_field: The magnetic field.
        isotropic_diffusion: Is the particle diffusion isotropic?
    """

    def __init__(self, magnetic_field, mfp, nr_steps, step_size):
        self.magnetic_field = magnetic_field
        self.mfp = mfp
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.isotropic_diffusion = False
        self.cartesian = False
        self.cylindrical = True

        self.init_jitclass_propagator()


    def get_description_propagator_type(self):
        print('propagation tpye: AnisotropicPropagator')