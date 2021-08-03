from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass

simulation_spec = [
    ('isotropic', b1),
    ('nr_steps', int32),
    ('dimensions', int32),
    ('step_distance', float32),
    ('step_size', float32),
    ('gyro_radius_eff', float32),
    ('speed', float32),

    ('chi_isotropic', float32),
    ('prob', float32[:]),

    ('pos', float32[:]),
    ('direction', float32[:]),
    ('phi', float32),
    ('distance', float32),
]

@jitclass(simulation_spec)
class Propagator():
    def __init__(self, nr_steps, step_size, mean_free_path):
        print('Propagator initialized')
        self.speed = 3*10**8 # [m^2/s]
        self.isotropic = False
        self.nr_steps = nr_steps
        self.step_size = step_size
        self.dimensions = 3
        self.pitch_angle_const = False
        
        xi = [self.speed / mean_free_path[0] / 2.0, self.speed / mean_free_path[1] / 2.0, self.speed / mean_free_path[2] / 2.0] # [1/s] frequency of change
        tau_step = self.step_size / self.speed
        self.prob = np.array([xi[0] * tau_step, xi[1] * tau_step, xi[2] * tau_step], dtype=np.float32)
        self.chi_isotropic = self.step_size / self.dimensions**0.5
        

    def change_direction(self, direction):
        ### change in direction happens with a propability that is defined by the 
        ### mean free path
        for p in range(self.dimensions):
            if np.random.random() < self.prob[p]:
                direction[p] = -1*direction[p]
        return direction


    def move_substep(self, pos, direction, phi, pitch_angle, distance, gyro_radius, s):
        distance = distance + self.step_size / self.dimensions
        self.gyro_radius_eff = gyro_radius / 3**0.5 # correcting for moving in rho direction (perp to phi) --> gyration increases by 2**0.5, which is why we have to divide here.
        if self.isotropic:
            pos = self.move_isotropic(pos, direction, pitch_angle, s)
        else:
            if s == 0:
                pos, phi = self.move_phi(pos, direction, phi, pitch_angle)
            if s == 1:
                pos, phi = self.move_rho(pos, direction, phi, pitch_angle)
            if s == 2:
                pos = self.move_isotropic(pos, direction, pitch_angle, 2)
        data = {
            'distance': distance, 
            'phi': phi,
            'pos': self.position(pos)
        }
        return data
         
            
    def move_isotropic(self, pos, direction, pitch_angle, s):
        if s == 2:
            distance_s = self.step_size * np.cos(pitch_angle)
        else:
            distance_s = self.step_size * np.sin(pitch_angle) / 2**0.5
        pos[s] = pos[s] + direction[s] * distance_s
        return self.position(pos)
        
        
    def move_phi(self, pos, direction, phi, pitch_angle):
        phi_old = phi
        phi = phi
        delta_rho = self.step_size * np.sin(pitch_angle)
        delta_phi = 2 * np.arcsin(delta_rho / (2 * 2**0.5 * self.gyro_radius_eff))
        phi = phi_old + delta_phi * direction[0]
        chi_x_1 = self.gyro_radius_eff * (np.cos(phi) - np.cos(phi_old))
        chi_y_1 = self.gyro_radius_eff * (np.sin(phi) - np.sin(phi_old))
        pos[0] = pos[0] + chi_x_1
        pos[1] = pos[1] + chi_y_1
        return self.position(pos), phi

                      
    def move_rho(self, pos, direction, phi, pitch_angle):
        delta_rho = self.step_size * np.sin(pitch_angle)
        chi_x_2 = np.cos(phi) * direction[1] * delta_rho
        chi_y_2 = np.sin(phi) * direction[1] * delta_rho
        pos[0] = pos[0] + chi_x_2
        pos[1] = pos[1] + chi_y_2
        return self.position(pos), phi


    def position(self, pos):
        return np.array([pos[0], pos[1], pos[2]], dtype=np.float32)