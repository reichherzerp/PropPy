"""The Particle State keeps track of all relevant particle properties.

In each simulation step, the particle state gets updated. It can be passed
to observers to easily write out the current particle state for analysis.
"""


from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass



particle_state_spec = [
    ('step_distance', float32),
    ('chi_isotropic', float32),
    ('speed', float32),
    ('energy', float32),
    ('isotropic_diffusion', b1),
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
    ('rad_prev', float32),
    ('direction', float32[:]),
    ('step', int32),
    ('substep', int32),
]

@jitclass(particle_state_spec)
class ParticleState():
    """Class for keeping track of the particle state.
    
    All properties of the particle have to be defined here. They will be initialized
    in the source and updated during the propagation process. In the observer module
    the particle state can be written out.

    Attributes:
        speed: Speed of the particle in [m^2/s].
        energy: Energy of the particle in [eV].
        gyroradius: Gyroradius of the particle in [m].
        gyroradius_eff: In the cylindircal movement, an effective radius is established that slightly differs from the above.
        particle_id: Unique id for easy analyzing them later.
        isotropic_diffusion: Parameter to determine if the diffusion is isotropic.
        dimensions: Number of dimensions.
        distance: Travelled distance of the particle.
        substep: Number of current substep.
        step: Number of current step. 
        pos_start: Initial position of the particle.
        pos: Current position of the particle.
        pos_prev: Previous position of the particle.
        rad_prev: Previous radius of the particle for spherical observer.
        direction: Directions for the random walk (can be 1 or -1 in each direction).
        dir: Direction into which particle points.
        phi: Angle in the xy-plane with respect to the x-axis.
        pitch_angle: Pitch angle of the particle.
    """

    def __init__(self, particle_id, energy, pos, phi, pitch_angle, dimensions):
        self.speed = 3*10**8 # [m^2/s]
        self.energy = energy
        self.gyroradius = 0.0
        self.gyroradius_eff = 0.0
        self.particle_id = particle_id
        self.isotropic_diffusion = False
        self.dimensions = dimensions
        self.distance = 0.0
        self.substep = 0
        self.step = 0
        self.pos_start = pos[:]
        self.pos = pos[:]
        self.pos_prev = self.pos[:]
        self.rad_prev = 0.0
        self.direction = np.array([1.0, 1.0, 1.0], dtype=np.float32)
        self.dir = np.array([np.cos(phi)*np.sin(pitch_angle), np.sin(phi)*np.sin(pitch_angle), np.cos(pitch_angle)], dtype=np.float32)
        self.phi = phi
        self.pitch_angle = pitch_angle


    def init_position(self):
        self.pos = np.array([self.pos_start[0], self.pos_start[1], self.pos_start[2]], dtype=np.float32)
        

