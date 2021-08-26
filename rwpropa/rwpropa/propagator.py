from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from .magnetic_field import *
from .particle_state import *
from abc import ABC, ABCMeta, abstractmethod


simulation_spec = [
    ('cartesian', b1),
    ('cylindrical', b1),
    ('nr_steps', int32),
    ('dimensions', int32),
    ('background_direction', int32),
    ('step_distance', float32),
    ('step_size', float32),
    ('gyroradius_eff', float32),
    ('speed', float32),
    ('pitch_angle_const', b1),

    ('prob', float32[:]),

    ('pos', float32[:]),
    ('direction', float32[:]),
    ('phi', float32),
    ('pitch_angle', float32),
    ('distance', float32),

    ('magnetic_field', MagneticField.class_type.instance_type),
    ('particle_state', ParticleState.class_type.instance_type),
]

@jitclass(simulation_spec)
class Propagator():
    def __init__(self, nr_steps, step_size, prob, magnetic_field):
        print('Propagator initialized')
        self.speed = 2.998*10**8 # [m/s]
        self.cartesian = False
        self.cylindrical = True
        self.nr_steps = nr_steps
        self.step_size = step_size # [m]
        self.dimensions = 3
        self.pitch_angle_const = True
        self.background_direction = 2 # direction of a background field
        self.magnetic_field = magnetic_field
        self.prob = prob
        

    def change_direction(self, direction):
        # change in direction happens with a propability that is defined by the 
        # mean free path and stored in self.prob
        for p in range(self.dimensions):
            if np.random.random() < self.prob[p]:
                direction[p] = -1*direction[p]
        return direction


    def change_pitch_angle(self, pitch_angle):
        # changes in the pitch angle are caused by resonant scattering of particles at
        # fluctuations of the turbulence that satisfy the resonance scattering criterion.
        # these changes in pitch angle are approximated in Kulsrud & Pearce (1969, ApJ, 156, 445) 
        # and Reichherzer et al. (2020, MNRAS) as follows:
        # delta mu = b/B.
        # here, b is the rms field strength of the turbulence and B the magnetic field strength of
        # the ordered magnetic field lines. Only valid for weak turbulence levels b << B.
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
        pitch_angle = pitch_angle_1
        return pitch_angle


    def propagate(self, particle_state):
        # 1. step
        self.compute_bfield_angles()
        # 2. step 
        self.compute_rotation_matrix()
        # 3. step
        self.transform_into_local_frame()
        # 4. step
        self.compute_local_angles()
        # 5. step
        self.update_directions()
        # 6 .step
        self.update_pitch_angle()
        # 7. step
        self.compute_local_movements()

        for substep in range(particle_state.dimensions):
            particle_state.substep = substep
            particle_state = self.move_substep(particle_state)
        
        return particle_state


    def compute_bfield_angles(self):
        # 1. alignment B vector and z axis
        # - compute angle theta_B between z-axis and B-field
        # - compute angle phi of B-filed
        # TODO: implement
        pass


    def compute_rotation_matrix(self):
        # 2. find rotation matrix to project vector B onto the z-axis
        # use angles theta_B and phi from first step
        # TODO: implement
        pass


    def transform_into_local_frame(self):
        # 3. rotate particle state into local frame
        # TODO: implement
        pass


    def compute_local_angles(self):
        # 4. compute pitch angle and phi in the local frame
        # (note that the pitch angle should be the same in the local and the global frame)
        # - compute phi
        # - compute pithc angle
        # TODO: implement
        pass


    def update_directions(self):
        # 5. update the directions 
        # should be function change_directio(self) from above
        # TODO: implement
        pass


    def update_pitch_angle(self):
        # 6. update the directions 
        # should be function change_pitch_angle(self) from above
        # TODO: implement
        pass


    def compute_local_movements(self):
        # 7. compute movement in local frame, where b field aligns z-axis
        # in 3d, this should be thre substeps
        # TODO: implement
        pass


    def move_substep(self, particle_state):
        # adapt phi and pitch angle in case the b-field vector changed: issue 26
        
        # division into local and global step is slow... TODO: check for performance optimization
        # local step
        particle_state, move_local_array = self.move_local(particle_state)
        # global step
        particle_state = self.move_global(particle_state, move_local_array)
        return particle_state


    def move_local(self, particle_state):
        if self.cartesian:
            # cartesian coordinates -> move in x, y and z directions
            particle_state, move_local = self.move_cartesian(particle_state)
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
        # find the roation matrix to roate magnetic field of local frame (0,0,1)
        # into global frame (self.magnetic_field.direction) 
        # follow procedure described in:
        # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
        # python: https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space

        # TODO: the following two function calls are really slow...
        rotation_matrix = self.find_rotation_matrix()
        move_global = self.rotate(rotation_matrix, move_local)
        for s in range(self.dimensions):
            # TODO: make transformation from local to global frame
            self.magnetic_field.direction
            ps.pos[s] = ps.pos[s] + move_local[s]
        return ps


    def find_rotation_matrix(self):
        rotation_matrix = [[0,0,0],[0,0,0],[0,0,0]]
        return rotation_matrix


    def rotate(self, rotation_matrix, move_local):
        # using the rotation matrix to transform the local movement 
        # into the global frame
        move_global = move_local
        return move_global
    
                 
    def move_cartesian(self, particle_state):
        distance_s = 0.0
        if particle_state.substep == self.background_direction:
            distance_s = self.step_size * np.cos(particle_state.pitch_angle) * particle_state.direction[particle_state.substep]
        else:
            distance_s = self.step_size * np.sin(particle_state.pitch_angle) / 2**0.5
        #particle_state.pos[particle_state.substep] = particle_state.pos[particle_state.substep] + particle_state.direction[particle_state.substep] * distance_s
        particle_state.distance = particle_state.distance + distance_s
        move_local = [0,0,0]
        for s in range(self.dimensions):
            if s == particle_state.substep:
                move_local[s] = distance_s
        return particle_state, self.float_array(move_local)
        
        
    def move_phi(self, particle_state):
        phi_old = particle_state.phi
        distance_s = self.step_size * np.sin(particle_state.pitch_angle) / 2**0.5
        particle_state.distance = particle_state.distance + distance_s
        delta_phi = self.compute_delta_phi(particle_state)
        particle_state.phi = phi_old + delta_phi * particle_state.direction[0]
        chi_x_1 = particle_state.gyroradius_eff * (np.cos(particle_state.phi) - np.cos(phi_old))
        chi_y_1 = particle_state.gyroradius_eff * (np.sin(particle_state.phi) - np.sin(phi_old))
        #particle_state.pos[0] = particle_state.pos[0] + chi_x_1
        #particle_state.pos[1] = particle_state.pos[1] + chi_y_1
        return particle_state, self.float_array([chi_x_1, chi_y_1, 0])

                      
    def move_rho(self, particle_state):
        distance_s = self.step_size * np.sin(particle_state.pitch_angle) / 2**0.5
        particle_state.distance = particle_state.distance + distance_s
        delta_rho = self.step_size * np.sin(particle_state.pitch_angle) / 2**0.5
        chi_x_2 = np.cos(particle_state.phi) * particle_state.direction[1] * delta_rho
        chi_y_2 = np.sin(particle_state.phi) * particle_state.direction[1] * delta_rho
        #particle_state.pos[0] = particle_state.pos[0] + chi_x_2
        #particle_state.pos[1] = particle_state.pos[1] + chi_y_2
        return particle_state, self.float_array([chi_x_2, chi_y_2, 0])


    def compute_delta_phi(self, ps):
        delta_rho = self.step_size * np.sin(ps.pitch_angle)
        delta_phi = 2 * np.arcsin(delta_rho / (2 * 2**0.5 * ps.gyroradius_eff))
        return delta_phi


    def float_array(self, pos):
        return np.array([pos[0], pos[1], pos[2]], dtype=np.float32)


    def set_gyroradius(self, ps):
        # default gyroradius for protons (v=c) with 1eV in magnetic field with strength 1Gaus
        gyroradius_0 = 3.336*10**(-5) # meters
        ps.gyroradius = gyroradius_0 * ps.energy / self.magnetic_field.rms # meters
        ps.gyroradius_eff = ps.gyroradius / 3**0.5 # correcting for moving in rho direction (perp to phi) --> gyration increases by 2**0.5, which is why we have to divide here. 
        return ps
        

    def set_pitch_angle_const(self, const_bool):
        # keep the pitch angle either constant or allow for changes 
        # during each propagation step.
        self.pitch_angle_const = const_bool


    def set_dimensions(self, dimensions):
        # default is 3d -> dimensions = 3
        # more than 3 dimensions are not supported
        self.dimensions = dimensions


    def set_cartesian_coords(self, cartesian):
        # there are cartesian or cylindrical coordinates available. 
        # cylindrical coordinates activated lead to the usage of
        # move_phi, move_rho and move_cartesian for the z-direction
        self.cartesian = cartesian
        self.cylindrical = not cartesian

    
    def set_cylindrical_coords(self, cylindrical):
        # there are cartesian or cylindrical coordinates available. 
        # cylindrical coordinates activated lead to the usage of
        # move_phi, move_rho and move_cartesian for the z-direction
        self.cartesian = not cylindrical
        self.cylindrical = cylindrical


    def set_speed(self, speed):
        # units = [m/s]
        # change the speed of the particles.
        # the default speed is the speed of light that is valid for
        # relativistic particles
        self.speed = speed


    def set_nr_steps(self, nr_steps):
        # change number of steps
        self.nr_steps = nr_steps

    
    def set_step_size(self, step_size):
        # units = [m]
        # change distance of each step that particles travel 
        self.nr_steps = step_size

    
    def set_prob_init(self, mean_free_path, speed, step_size):
        xi = [
            speed / mean_free_path[0] / 2.0, 
            speed / mean_free_path[1] / 2.0, 
            speed / mean_free_path[2] / 2.0
        ] # [1/s] frequency of change
        tau_step = step_size / speed
        return np.array([xi[0] * tau_step, xi[1] * tau_step, xi[2] * tau_step], dtype=np.float32)


    def set_prob(self, mean_free_path):
        self.prob = self.set_prob_init(mean_free_path, self.propagator.speed, self.propagator.step_size)
    

    def set_magnetic_field(self, magnetic_field):
        self.magnetic_field = magnetic_field


    def get_description(self):
        # note: description does not contain the information of the underling special propagator 
        # (if there was one)
        # that was used during adding the propagator to the simulation. 
        # To get this info, get_description(self) 
        # has to be called directly on the instance of the used propagator (see tutorials for details)
        self.get_description_general()
        self.get_description_parameters()


    def get_description_general(self):
        # print out the discription of the object with all relevant instance parameters
        print("""Description Propagator:
                The propagator is responsible for the movement of the particles. 
                It performs the change of direction and the movement in the respective direction.
                There are two phases:
                 - change direction with probability (see below)
                 - move in all directions
                The movement takes place according to the random walk (RW).\n""")


    def get_description_parameters(self):
        # print out the discription of the object with all relevant instance parameters
        print("""Description Propagator:
                The propagator is responsible for the movement of the particles. 
                It performs the change of direction and the movement in the respective direction.
                There are two phases:
                 - change direction with probability (see below)
                 - move in all directions
                The movement takes place according to the random walk (RW).\n""")
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
        print('particle speed: ' ,self.speed, ' m/s')
        print('number steps: ', self.nr_steps)  
        print('step size: ', self.step_size, ' m')  
        print('step duration: ', self.step_size / self.speed, ' s') 
        print('total distance: ', self.step_size * self.nr_steps, ' m')
        print('total duration: ', self.step_size * self.nr_steps / self.speed, ' s')
        print('probability to change directions in step: ', self.prob*100, '%')  




#---------------------------------------------------------------------------------
# below are the abstract base class and all sub classes of the special propagators
# that have to be added to the simulation. Each special propagator stores a
# Propagator object in its instance parameter to be used in the simulation
# This diversions is needed because numba does not support 
# inheritance via ABC and Propagator() needs the label @jitclass as it is called 
# during the numba optimized simulation loop of the run_simulation() function. 
# This workaround supports both concepts with the 
# advantages of fast code and easy addition of new propagator where the structure 
# is now defined by the Abstract Base class and enforeced via the ABCMeta class



class AbstractPropagatorMeta(ABCMeta):
    # required attributes that have to be implemented in __init__ of all
    # sub classes
    required_attributes = []

    def __call__(self, *args, **kwargs):
        # check if required attributes that have to be implemented in __init__ of all
        # sub classes are really implemented. Raise an error if not
        obj = super(AbstractPropagatorMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj



class AbstractPropagator(object, metaclass=AbstractPropagatorMeta):
    # abstract base class for all special observers.
    # functions with the label @abstractmethod have to be implemented in 
    # the special observer classes

    # all required_attributes have to be implemented in sub classes
    required_attributes = [
        'propagator', 
        'dimensions', 
        'speed', 
        'mfp', 
        'nr_steps', 
        'magnetic_field', 
        'step_size', 
        'isotropic'
    ]
 
    @abstractmethod
    def __init__(self, order):
        # implementation required in all sub classes.
        # all required_attributes have to be implemented in sub classes
        pass

    def set_basic_parameters(self):
        self.dimensions = 3
        self.speed = 2.998*10**8 # [m/s]


    def set_pitch_angle_const(self, const_bool):
        # keep the pitch angle either constant or allow for changes 
        # during each propagation step.
        self.pitch_angle_const = const_bool
        self.propagator.pitch_angle_const = const_bool


    def set_dimensions(self, dimensions):
        # default is 3d -> dimensions = 3
        # more than 3 dimensions are not supported
        self.dimensions = dimensions
        self.propagator.dimensions = dimensions


    def set_cartesian_coords(self, cartesian):
        # there are cartesian or cylindrical coordinates available. 
        # cylindrical coordinates activated lead to the usage of
        # move_phi, move_rho and move_cartesian for the z-direction
        self.cartesian = cartesian
        self.propagator.cylindrical = not cartesian

    
    def set_cylindrical_coords(self, cylindrical):
        # there are cartesian or cylindrical coordinates available. 
        # cylindrical coordinates activated lead to the usage of
        # move_phi, move_rho and move_cartesian for the z-direction
        self.cartesian = not cylindrical
        self.cylindrical = cylindrical
        self.propagator.cartesian = not cylindrical
        self.propagator.cylindrical = cylindrical


    def set_speed(self, speed):
        # units = [m/s]
        # change the speed of the particles.
        # the default speed is the speed of light that is valid for
        # relativistic particles
        self.speed = speed
        self.propagator.speed = speed


    def set_nr_steps(self, nr_steps):
        # change number of steps
        self.nr_steps = nr_steps
        self.propagator.nr_steps = nr_steps

    
    def set_step_size(self, step_size):
        # units = [m]
        # change distance of each step that particles travel 
        self.nr_steps = step_size
        self.propagator.nr_steps = step_size

    
    def set_prob_init(self, mean_free_path, speed, step_size):
        xi = [speed / mean_free_path[0] / 2.0, speed / mean_free_path[1] / 2.0, speed / mean_free_path[2] / 2.0] # [1/s] frequency of change
        tau_step = step_size / speed
        return np.array([xi[0] * tau_step, xi[1] * tau_step, xi[2] * tau_step], dtype=np.float32)


    def set_prob(self, mean_free_path):
        self.propagator.prob = self.set_prob_init(mean_free_path, self.propagator.speed, self.propagator.step_size)
        self.prob = self.propagator.prob
    

    def set_magnetic_field(self, magnetic_field):
        # allow to set the propagator magnetic field
        self.magnetic_field = magnetic_field
        self.propagator.magnetic_field = magnetic_field


    def convert_mfp_input(self, mfp_input):
        # check the input of the mean free paths
        if self.isotropic == False:
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

        return mfp

    
    def init_jitclass_propagator(self):
        self.set_basic_parameters()
        self.mfp = self.convert_mfp_input(self.mfp)
        mfp_final = self.set_prob_init(self.mfp, self.speed, self.step_size)
        propagator = Propagator(self.nr_steps, self.step_size, mfp_final, self.magnetic_field)
        self.propagator = propagator

    
    @abstractmethod
    def get_description_propagator_type(sefl):
        pass


    def get_description(self):
        # print the information of the relevant parameters and the description of 
        # the special propagation type that was chosen
        self.propagator.get_description_general()
        self.get_description_propagator_type()
        self.propagator.get_description_parameters()


    def get_description_general(self):
        # called by all special propagator classes below.
        # introduction of the description output
        print("""Description Propagator:
                The propagator is responsible for the movement of the particles. 
                It performs the change of direction and the movement in the respective direction.
                There are two phases:
                 - change direction with probability (see below)
                 - move in all directions
                The movement takes place according to the random walk (RW).\n""")


    def get_description_parameters(self):   
        # called by all special propagator classes below.
        # print out all relevant instance parameters
        print('particle speed: ' ,self.speed, ' m/s')
        print('number steps: ', self.nr_steps)  
        print('step size: ', self.step_size, ' m')  
        print('step duration: ', self.step_size / self.speed, ' s') 
        print('total distance: ', self.step_size * self.nr_steps, ' m')
        print('total duration: ', self.step_size * self.nr_steps / self.speed, ' s') 
        print('call get_description directly on the propagator that was added to the simulation:\n')
        print('sim = rwpropa.Simulation()')
        print('...')
        print('sim.add_propagator(propagator)')
        print('sim.propagator.get_description()')



class IsotropicPropagatorDefault(AbstractPropagator):
    def __init__(self):
        self.nr_steps = 2*10**5
        self.step_size = 0.5*10**10 # [m]
        # isotropic diffusion coefficient
        self.mfp = np.array([10**12, 10**12, 10**12], dtype=np.float32)  # [m]
        # no background magnetic field
        self.magnetic_field = OrderedBackgroundField(0, [0,0,1]).magnetic_field
        self.isotropic = True

        self.init_jitclass_propagator() 


    def get_description_propagator_type(self):
        print('propagation tpye: IsotropicPropagatorDefault')



class IsotropicPropagator(AbstractPropagator):
    def __init__(self, magnetic_field, mfp, nr_steps, step_size):
        self.magnetic_field = magnetic_field
        self.mfp = mfp
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.isotropic = True

        self.init_jitclass_propagator()


    def get_description_propagator_type(self):
        print('propagation tpye: IsotropicPropagator')
  


class AnisotropicPropagator(AbstractPropagator):
    def __init__(self, magnetic_field, mfp, nr_steps, step_size):
        self.magnetic_field = magnetic_field
        self.mfp = mfp
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.isotropic = False

        self.init_jitclass_propagator()


    def get_description_propagator_type(self):
        print('propagation tpye: AnisotropicPropagator')