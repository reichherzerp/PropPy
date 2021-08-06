from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from .magnetic_field import *
from abc import ABC, ABCMeta, abstractmethod


simulation_spec = [
    ('cartesian', b1),
    ('cylindrical', b1),
    ('nr_steps', int32),
    ('dimensions', int32),
    ('background_direction', int32),
    ('step_distance', float32),
    ('step_size', float32),
    ('gyro_radius_eff', float32),
    ('speed', float32),
    ('pitch_angle_const', b1),

    ('prob', float32[:]),

    ('pos', float32[:]),
    ('direction', float32[:]),
    ('phi', float32),
    ('pitch_angle', float32),
    ('distance', float32),

    ('magnetic_field', MagneticField.class_type.instance_type),
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


    def move_substep(self, pos, direction, phi, pitch_angle, distance, gyro_radius, s):
        self.gyro_radius_eff = gyro_radius / 3**0.5 # correcting for moving in rho direction (perp to phi) --> gyration increases by 2**0.5, which is why we have to divide here.
        if self.cartesian:
            # cartesian coordinates -> move in x, y and z directions
            pos, distance = self.move_cartesian(pos, direction, pitch_angle, distance, s)
        else:
            # cylindrical coordinates -> move in phi, rho and z directions
            if s == 0:
                pos, phi, distance = self.move_phi(pos, direction, phi, pitch_angle, distance)
            if s == 1:
                pos, phi, distance = self.move_rho(pos, direction, phi, pitch_angle, distance)
            if s == 2:
                pos, distance = self.move_cartesian(pos, direction, pitch_angle, distance, 2)
        data = {
            'distance': distance, 
            'phi': phi,
            'pos': self.position(pos)
        }
        return data
         
            
    def move_cartesian(self, pos, direction, pitch_angle, distance, s):
        if s == self.background_direction:
            distance_s = self.step_size * np.cos(pitch_angle)
            if self.pitch_angle_const == False:
                pos[s] = pos[s] + distance_s
            else:
                pos[s] = pos[s] + direction[s] * distance_s
        else:
            distance_s = self.step_size * np.sin(pitch_angle) / 2**0.5
            pos[s] = pos[s] + direction[s] * distance_s
        distance = distance + distance_s
        return self.position(pos), distance
        
        
    def move_phi(self, pos, direction, phi, pitch_angle, distance):
        phi_old = phi
        distance_s = self.step_size * np.sin(pitch_angle) / 2**0.5
        distance = distance + distance_s
        delta_phi = self.compute_delta_phi(pitch_angle)
        phi = phi_old + delta_phi * direction[0]
        chi_x_1 = self.gyro_radius_eff * (np.cos(phi) - np.cos(phi_old))
        chi_y_1 = self.gyro_radius_eff * (np.sin(phi) - np.sin(phi_old))
        pos[0] = pos[0] + chi_x_1
        pos[1] = pos[1] + chi_y_1
        return self.position(pos), phi, distance

                      
    def move_rho(self, pos, direction, phi, pitch_angle, distance):
        distance_s = self.step_size * np.sin(pitch_angle) / 2**0.5
        distance = distance + distance_s
        delta_rho = self.step_size * np.sin(pitch_angle) / 2**0.5
        chi_x_2 = np.cos(phi) * direction[1] * delta_rho
        chi_y_2 = np.sin(phi) * direction[1] * delta_rho
        pos[0] = pos[0] + chi_x_2
        pos[1] = pos[1] + chi_y_2
        return self.position(pos), phi, distance


    def compute_delta_phi(self, pitch_angle):
        delta_rho = self.step_size * np.sin(pitch_angle)
        delta_phi = 2 * np.arcsin(delta_rho / (2 * 2**0.5 * self.gyro_radius_eff))
        return delta_phi


    def position(self, pos):
        return np.array([pos[0], pos[1], pos[2]], dtype=np.float32)


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
        xi = [speed / mean_free_path[0] / 2.0, speed / mean_free_path[1] / 2.0, speed / mean_free_path[2] / 2.0] # [1/s] frequency of change
        tau_step = step_size / speed
        return np.array([xi[0] * tau_step, xi[1] * tau_step, xi[2] * tau_step], dtype=np.float32)


    def set_prob(self, mean_free_path):
        self.prob = self.set_prob_init(mean_free_path, self.propagator.speed, self.propagator.step_size)
    

    def set_magnetic_field(self, magnetic_field):
        self.magnetic_field = magnetic_field


    def get_description(self):
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




#--------------------------------------------------------------------------------------------------------------






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



class IsotropicPropagator(AbstractPropagator):
    def __init__(self, magnetic_field, mfp, nr_steps, step_size):
        self.magnetic_field = magnetic_field
        self.mfp = mfp
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.isotropic = True

        self.init_jitclass_propagator()
  



class AnisotropicPropagator(AbstractPropagator):
    def __init__(self, magnetic_field, mfp, nr_steps, step_size):
        self.magnetic_field = magnetic_field
        self.mfp = mfp
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.isotropic = False

        self.init_jitclass_propagator()














