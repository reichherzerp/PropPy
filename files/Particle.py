from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from files.Observer import Observer

simulation_spec = [
    ('step_distance', float32),
    ('chi_isotropic', float32),
    ('speed', float32),
    ('isotropic', b1),
    ('distance', float32),
    ('gyro_radius_eff', float32),
    ('phi', float32),
    ('particle_id', int32),
    ('dimensions', int32),
    ('pos_start', float32[:]),
    ('pos', float32[:]),
    ('pos_prev', float32[:]),
    ('direction', float32[:]),
    ('prob', float32[:]),
    ('observer', Observer.class_type.instance_type),
]

@jitclass(simulation_spec)
class Particle():
    def __init__(self, particle_id, gyro_radius, mean_free_path, pos):
        self.speed = 3*10**8 # [m^2/s]
        self.step_distance = 0.5*10**10 ## [m]
        self.gyro_radius_eff = gyro_radius / 3**0.5 # correcting for moving in rho direction (perp to phi) --> gyration increases by 2**0.5, which is why we have to divide here.
        self.particle_id = particle_id
        self.isotropic = False
        self.dimensions = 3
        self.distance = 0.0
        self.pos_start = pos[:]
        self.pos = pos[:]
        self.pos_prev = self.pos[:]
        self.direction = np.array([1.0, 1.0, 1.0], dtype=np.float32)
        self.phi = 0.0
        xi = [self.speed / mean_free_path[0] / 2.0, self.speed / mean_free_path[1] / 2.0, self.speed / mean_free_path[2] / 2.0] # [1/s] frequency of change
        tau_step = self.step_distance / self.speed
        self.prob = np.array([xi[0] * tau_step, xi[1] * tau_step, xi[2] * tau_step], dtype=np.float32)
        self.chi_isotropic = self.step_distance / self.dimensions**0.5
        
        
    def simulate(self, observer, nr_steps):
        simulation_data = []
    
        simulation_data.append([self.particle_id, 0, self.distance, self.pos[0], self.pos[1], self.pos[2], -1.0])
        self.pos = np.array([self.pos_start[0], self.pos_start[1], self.pos_start[2]], dtype=np.float32)
        for i in range(1, nr_steps): 
            self.change_direction()
            self.pos_prev = self.pos 
            for s in range(self.dimensions):
                self.move_substep(s)
                observation = observer.observe(i, s, self.distance, self.pos, self.particle_id)
                if observation is not None:
                    simulation_data.append(observation)
                
        return simulation_data
                
        
    def change_direction(self):
        ### change in direction happens with a propability that is defined by the 
        ### mean free path
        for p in range(self.dimensions):
            if np.random.random() < self.prob[p]:
                self.direction[p] = -1*self.direction[p]
    

    def move_substep(self, s):
        self.distance = self.distance + self.step_distance / self.dimensions
        if self.isotropic:
            self.move_isotropic(s)
        else:
            if s == 0:
                self.move_phi()
            if s == 1:
                self.move_rho()
            if s == 2:
                self.move_isotropic(2)
            
            
    def move_isotropic(self, s):
        self.pos[s] = self.pos[s] + self.direction[s] * self.chi_isotropic
        
        
    def move_phi(self):
        distance_in_step = self.step_distance
        phi_old = self.phi
        phi = self.phi
        delta_phi = 2 * np.arcsin(distance_in_step / (12**0.5 * self.gyro_radius_eff))
        phi = phi_old + delta_phi * self.direction[0]
        self.phi = phi
        chi_x_1 = self.gyro_radius_eff * (np.cos(phi) - np.cos(phi_old)) 
        chi_y_1 = self.gyro_radius_eff * (np.sin(phi) - np.sin(phi_old))
        self.pos[0] = self.pos[0] + chi_x_1
        self.pos[1] = self.pos[1] + chi_y_1

                      
    def move_rho(self):
        distance_in_step = self.step_distance
        phi = self.phi
        delta_rho = distance_in_step / 3**0.5
        chi_x_2 = np.cos(phi) * self.direction[1] * delta_rho
        chi_y_2 = np.sin(phi) * self.direction[1] * delta_rho

        self.pos[0] = self.pos[0] + chi_x_2
        self.pos[1] = self.pos[1] + chi_y_2