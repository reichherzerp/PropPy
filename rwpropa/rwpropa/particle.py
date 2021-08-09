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
    ('particle_state', ParticleState.class_type.instance_type),
]

@jitclass(particle_spec)
class Particle():
    def __init__(self, particle_id, gyro_radius, pos, phi, pitch_angle, dimensions):
        self.particle_state = ParticleState(particle_id, gyro_radius, pos, phi, pitch_angle, dimensions)
        
        
    def simulate(self, observer, propagator):
        simulation_data = []
    
        simulation_data.append([self.particle_state.particle_id, 0, self.particle_state.distance, self.particle_state.pos[0], self.particle_state.pos[1], self.particle_state.pos[2], -1.0, self.particle_state.dimensions-1])
        self.particle_state.init_position()
        for i in range(1, propagator.nr_steps): 
            self.particle_state.direction = propagator.change_direction(self.particle_state.direction)
            self.particle_state.pitch_angle = propagator.change_pitch_angle(self.particle_state.pitch_angle)
            self.particle_state.pos_prev = self.particle_state.pos 
            for substep in range(self.particle_state.dimensions):
                self.particle_state.substep = substep
                self.propagate(propagator)
                observation = observer.observe(i, substep, self.particle_state.distance, self.particle_state.pos, self.particle_state.particle_id, self.particle_state.phi, self.particle_state.pitch_angle)
                if observation is not None:
                    simulation_data.append(observation)
                
        return simulation_data


    def propagate(self, propagator):
        self.particle_state = propagator.move_substep(self.particle_state)



