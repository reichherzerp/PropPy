from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass



particle_state_spec = [
    ('step_distance', float32),
    ('chi_isotropic', float32),
    ('speed', float32),
    ('energy', float32),
    ('isotropic', b1),
    ('distance', float32),
    ('gyroradius', float32),
    ('gyroradius_eff', float32),
    ('local_move', float32[:]),
    ('phi', float32),
    ('pitch_angle', float32),
    ('particle_id', int32),
    ('dimensions', int32),
    ('pos_start', float32[:]),
    ('pos', float32[:]),
    ('dir', float32[:]),
    ('pos_prev', float32[:]),
    ('direction', float32[:]),
    ('step', int32),
    ('substep', int32),
]

@jitclass(particle_state_spec)
class ParticleState():
    def __init__(self, particle_id, energy, pos, phi, pitch_angle, dimensions):
        self.speed = 3*10**8 # [m^2/s]
        self.energy = energy
        self.gyroradius = 0.0
        self.gyroradius_eff = 0.0
        self.particle_id = particle_id
        self.isotropic = False
        self.dimensions = dimensions
        self.distance = 0.0
        self.substep = 0
        self.step = 0
        self.pos_start = pos[:]
        self.pos = pos[:]
        self.pos_prev = self.pos[:]
        # directions for the random walk
        self.direction = np.array([1.0, 1.0, 1.0], dtype=np.float32)
        # direction into which particle points
        self.dir = np.array([np.cos(phi)*np.sin(pitch_angle), np.sin(phi)*np.sin(pitch_angle), np.cos(pitch_angle)], dtype=np.float32)
        self.phi = phi
        self.pitch_angle = pitch_angle


    def init_position(self):
        self.pos = np.array([self.pos_start[0], self.pos_start[1], self.pos_start[2]], dtype=np.float32)
        

