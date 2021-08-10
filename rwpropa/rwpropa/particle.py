from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from .observer import *
from .propagator import *
from .particle_state import *



particle_spec = [
    ('step_distance', float32),
    ('chi_isotropic', float32),
    ('speed', float32),
    ('isotropic', b1),
    ('distance', float32),
    ('gyro_radius', float32),
    ('phi', float32),
    ('pitch_angle', float32),
    ('particle_id', int32),
    ('dimensions', int32),
    ('pos_start', float32[:]),
    ('pos', float32[:]),
    ('pos_prev', float32[:]),
    ('direction', float32[:]),
    ('prob', float32[:]),
    ('observer', Observer.class_type.instance_type),
    ('propagator', Propagator.class_type.instance_type),
    ('ps', ParticleState.class_type.instance_type),
]

@jitclass(particle_spec)
class Particle():
    def __init__(self, particle_id, gyro_radius, pos, phi, pitch_angle, dimensions):
        self.ps = ParticleState(particle_id, gyro_radius, pos, phi, pitch_angle, dimensions)
        
        
    def simulate(self, observer, propagator):
        simulation_data = []
        simulation_data.append(observer.data_row(0, self.ps))
        self.ps.init_position()
        for i in range(1, propagator.nr_steps): 
            self.ps.step = i
            self.ps.pos_prev = self.ps.pos
            self.ps.direction = propagator.change_direction(self.ps.direction)
            self.ps.pitch_angle = propagator.change_pitch_angle(self.ps.pitch_angle) 
            for substep in range(self.ps.dimensions):
                self.ps.substep = substep
                self.propagate(propagator)
                observation = observer.observe(i, self.ps)
                if observation is not None:
                    simulation_data.append(observation)
                
        return simulation_data


    def propagate(self, propagator):
        self.ps = propagator.move_substep(self.ps)


    





