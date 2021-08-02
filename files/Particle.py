from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from files.Observer import Observer
from files.Propagator import Propagator

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
    ('propagator', Propagator.class_type.instance_type),
]

@jitclass(simulation_spec)
class ParticleState():
    def __init__(self):
        self.pos = np.array([1.0, 1.0, 1.0], dtype=np.float32)
        self.direction = np.array([1.0, 1.0, 1.0], dtype=np.float32)
        self.phi = 0.0


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
    ('propagator', Propagator.class_type.instance_type),
]

@jitclass(simulation_spec)
class Particle():
    def __init__(self, particle_id, gyro_radius, pos):
        self.speed = 3*10**8 # [m^2/s]
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

        
        
    def simulate(self, observer, propagator):
        simulation_data = []
    
        simulation_data.append([self.particle_id, 0, self.distance, self.pos[0], self.pos[1], self.pos[2], -1.0, self.dimensions-1])
        self.pos = np.array([self.pos_start[0], self.pos_start[1], self.pos_start[2]], dtype=np.float32)
        for i in range(1, propagator.nr_steps): 
            self.direction = propagator.change_direction(self.direction)
            self.pos_prev = self.pos 
            for substep in range(self.dimensions):
                self.distance, self.pos = propagator.move_substep(self.pos, self.direction, self.phi, self.distance, substep)
                
                observation = observer.observe(i, substep, self.distance, self.pos, self.particle_id)
                if observation is not None:
                    simulation_data.append(observation)
                
        return simulation_data
                
        
    
    
