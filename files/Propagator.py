from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass

simulation_spec = [
    ('cartesian', b1),
    ('cylindrical', b1),
    ('nr_steps', int32),
    ('dimensions', int32),
    ('step_distance', float32),
    ('step_size', float32),
    ('gyro_radius_eff', float32),
    ('speed', float32),
    ('pitch_angle_const', b1),

    ('chi_isotropic', float32),
    ('prob', float32[:]),

    ('pos', float32[:]),
    ('direction', float32[:]),
    ('phi', float32),
    ('pitch_angle', float32),
    ('distance', float32),
]

@jitclass(simulation_spec)
class Propagator():
    def __init__(self, nr_steps, step_size, mean_free_path):
        print('Propagator initialized')
        self.speed = 3*10**8 # [m^2/s]
        self.cartesian = False
        self.cylindrical = True
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.dimensions = 3
        self.pitch_angle_const = True
        
        xi = [self.speed / mean_free_path[0] / 2.0, self.speed / mean_free_path[1] / 2.0, self.speed / mean_free_path[2] / 2.0] # [1/s] frequency of change
        tau_step = self.step_size / self.speed
        self.prob = np.array([xi[0] * tau_step, xi[1] * tau_step, xi[2] * tau_step], dtype=np.float32)
        self.chi_isotropic = self.step_size / self.dimensions**0.5


    def set_pitch_angle_const(self, const_bool):
        self.pitch_angle_const = const_bool


    def set_dimensions(self, dimensions):
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
        self.speed = speed


    def change_direction(self, direction, pitch_angle):
        ### change in direction happens with a propability that is defined by the 
        ### mean free path
        for p in range(self.dimensions):
            if np.random.random() < self.prob[p]:
                if self.pitch_angle_const == False and p == self.dimensions-1:
                    if np.random.random() < 0.5:
                        direction[p] = -1*direction[p]
                    pitch_angle = pitch_angle + direction[2] * 0.1
                else:
                    direction[p] = -1*direction[p]
        return direction, pitch_angle


    def move_substep(self, pos, direction, phi, pitch_angle, distance, gyro_radius, s):
        distance = distance + self.step_size / self.dimensions
        self.gyro_radius_eff = gyro_radius / 3**0.5 # correcting for moving in rho direction (perp to phi) --> gyration increases by 2**0.5, which is why we have to divide here.
        if self.cartesian:
            # cartesian coordinates -> move in x, y and z directions
            pos = self.move_cartesian(pos, direction, pitch_angle, s)
        else:
            # cylindrical coordinates -> move in phi, rho and z directions
            if s == 0:
                pos, phi = self.move_phi(pos, direction, phi, pitch_angle)
            if s == 1:
                pos, phi = self.move_rho(pos, direction, phi, pitch_angle)
            if s == 2:
                pos = self.move_cartesian(pos, direction, pitch_angle, 2)
        data = {
            'distance': distance, 
            'phi': phi,
            'pitch_angle': pitch_angle,
            'pos': self.position(pos)
        }
        return data
         
            
    def move_cartesian(self, pos, direction, pitch_angle, s):
        if s == 2:
            distance_s = self.step_size * np.cos(pitch_angle)
            if self.pitch_angle_const == False:
                pos[s] = pos[s] + distance_s
            else:
                pos[s] = pos[s] + direction[s] * distance_s
        else:
            distance_s = self.step_size * np.sin(pitch_angle) / 2**0.5
            pos[s] = pos[s] + direction[s] * distance_s
        return self.position(pos)
        
        
    def move_phi(self, pos, direction, phi, pitch_angle):
        phi_old = phi
        delta_phi = self.compute_delta_phi(pitch_angle)
        phi = phi_old + delta_phi * direction[0]
        chi_x_1 = self.gyro_radius_eff * (np.cos(phi) - np.cos(phi_old))
        chi_y_1 = self.gyro_radius_eff * (np.sin(phi) - np.sin(phi_old))
        pos[0] = pos[0] + chi_x_1
        pos[1] = pos[1] + chi_y_1
        return self.position(pos), phi

                      
    def move_rho(self, pos, direction, phi, pitch_angle):
        delta_rho = self.step_size * np.sin(pitch_angle) / 2**0.5
        chi_x_2 = np.cos(phi) * direction[1] * delta_rho
        chi_y_2 = np.sin(phi) * direction[1] * delta_rho
        pos[0] = pos[0] + chi_x_2
        pos[1] = pos[1] + chi_y_2
        return self.position(pos), phi


    def compute_delta_phi(self, pitch_angle):
        delta_rho = self.step_size * np.sin(pitch_angle)
        delta_phi = 2 * np.arcsin(delta_rho / (2 * 2**0.5 * self.gyro_radius_eff))
        return delta_phi


    def position(self, pos):
        return np.array([pos[0], pos[1], pos[2]], dtype=np.float32)