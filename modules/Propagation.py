import random
from numba import jit, int32, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass
import numpy as np

spec = [   
    ('step_size', float32), # a simple scalar field
    ('chi_normalized', float32[:]),
    ('nr_fixed_gyro_radii', int32),
]

@jitclass(spec)
class Propagation():
    def __init__(self, step_size):
        self.step_size = step_size
        self.chi_normalized = np.array([1.0/3**0.5, 1.0/3**0.5, 1.0/3**0.5], dtype=np.float32)
        self.nr_fixed_gyro_radii = 1

    def set_nr_fixed_gyro_radii(self, nr_fixed_gyro_radii):
        self.nr_fixed_gyro_radii = nr_fixed_gyro_radii

    def set_normalized_chi(self, chi_normalized):
        self.chi_normalized = chi_normalized
        
    def move(self, particle):
        pos_previous = particle.pos
        kappa = particle.diffusion_tensor
        tau_step = 1.0
        
        v = 1.0
        chi = self.chi_normalized*v
        tau = kappa*tau_step**2/chi**2
        xi = 2/tau
        prop_turn = tau_step*xi 
        
        step_size = 1
        particle.t = particle.t + 1*step_size
        
        
        ### 1. change direction (if needed)
        random_value = random.random()
        if particle.t > self.nr_fixed_gyro_radii * particle.gyro_radius:
            # after gyroradius ~ half gyroorbit, particle 
            # needs to change direction perp to background field
            if (random_value <= prop_turn[0]):
                particle.direction[0] = particle.direction[0]*(-1)
            if (random_value <= prop_turn[1]):
                particle.direction[1] = particle.direction[1]*(-1)
        if (random_value <= prop_turn[2]):
            particle.direction[2] = particle.direction[2]*(-1)
          
        ### 2. move along direction
        abs_direction = (particle.direction[0]**2+particle.direction[1]**2+particle.direction[2]**2)**0.5
        speed = particle.direction / abs_direction * v
        d = speed * tau_step 
        particle.pos[0] = pos_previous[0] + d[0]
        particle.pos[1] = pos_previous[1] + d[1]
        particle.pos[2] = pos_previous[2] + d[2]
        particle.pos_previous[0] = pos_previous[0] - d[0]
        particle.pos_previous[1] = pos_previous[1] - d[1]
        particle.pos_previous[2] = pos_previous[2] - d[2]


    def move_polar(self, particle):
        pos_previous = particle.pos
        kappa = particle.diffusion_tensor
        tau_step = 1.0
        
        v = 1.0
        chi = self.chi_normalized*v
        tau = kappa*tau_step**2/chi**2
        xi = 2/tau
        prop_turn = tau_step*xi 
        
        step_size = 1
        particle.t = particle.t + 1*step_size
        
        
        ### 1. change direction (if needed)
        random_value = random.random()
        if particle.t > self.nr_fixed_gyro_radii * particle.gyro_radius:
            # after gyroradius ~ half gyroorbit, particle 
            # needs to change direction perp to background field
            if (random_value <= prop_turn[0]):
                particle.direction[0] = particle.direction[0]*(-1)
            if (random_value <= prop_turn[1]):
                particle.direction[1] = particle.direction[1]*(-1)
        if (random_value <= prop_turn[2]):
            particle.direction[2] = particle.direction[2]*(-1)
          
        ### 2. move along direction
        abs_direction = (particle.direction[0]**2+particle.direction[1]**2+particle.direction[2]**2)**0.5
        speed = particle.direction / abs_direction * v
        d = speed * tau_step 
        ### move along phi
        particle.phi = particle.phi + particle.direction[1] * 0.05 / 2 / np.pi
        ### move along r:
        delta_x = particle.gyro_radius * np.cos(particle.phi) * particle.direction[0]
        delta_y = particle.gyro_radius * np.sin(particle.phi) * particle.direction[0]
        r2 = (delta_x**2+delta_y**2)**0.5
        particle.pos[0] = pos_previous[0] + delta_x/r2
        particle.pos[1] = pos_previous[1] + delta_y/r2
        particle.pos[2] = pos_previous[2] + particle.direction[2]
        particle.pos_previous[0] = pos_previous[0] - delta_x/r2
        particle.pos_previous[1] = pos_previous[1] - delta_y/r2
        particle.pos_previous[2] = pos_previous[2] - particle.direction[2]
 
