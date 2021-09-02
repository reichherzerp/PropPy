"""Individual particles will be simulated here.

By calling the particles simulate function, the particle will be
propagated completely and observed as specifyed in the observer.
"""

from numba import jit, b1, float32, int32
import numpy as np
from numba.experimental import jitclass
from .observer import *
from .propagator import *
from .particle_state import *



particle_spec = [
    ('observer', Observer.class_type.instance_type),
    ('propagator', Propagator.class_type.instance_type),
    ('ps', ParticleState.class_type.instance_type),
]

@jitclass(particle_spec)
class Particle():
    """Particle class for simulation particles. 
     
    The properties of the particles are stored within the particle state object.
    By calling the particles simulate function, the particle will be
    propagated completely and observed as specifyed in the observer.

    Attributes:
        observer: User specifyed observer.
        propagator: User specifyed propagator.
        ps: ParticleState of the current particle.
 
    """
    def __init__(self, particle_id, energy, pos, phi, pitch_angle, dimensions):
        self.ps = ParticleState(particle_id, energy, pos, phi, pitch_angle, dimensions)
        
        
    def simulate(self, observer, propagator):
        """Core simulation function of a single particle.
        
        The particle is propagated as specifyed in the list of steps in the propagator.
        In each propagation step a particle may be observed, depending on the observation
        conditions.

        Args:
            observer: User specifyed observer.
            propagator: User specifyed propagator.
        Returns:
            simulation_data: Observed data of the particles.
        """
        simulation_data = []

        self.ps.init_position()
        for step in range(0, propagator.nr_steps+1): 
            self.start(propagator, step)
            self.ps.pos_prev =  np.array([self.ps.pos[0], self.ps.pos[1], self.ps.pos[2]], dtype=np.float32)
            for substep in range(self.ps.dimensions):
                self.propagate_substep(propagator, substep)
                observation = observer.observe(self.ps)
                if observation is not None:
                    simulation_data.append(observation)

                    
        return simulation_data


    def start(self, propagator, step):
        """Start the simulation by using the parameters defined in the source.
        """
        self.ps.step = step
        self.ps.pos_prev = self.ps.pos
        self.ps = propagator.set_gyroradius(self.ps)
        self.ps.direction = propagator.change_direction(self.ps.direction)
        self.ps.pitch_angle = propagator.change_pitch_angle(self.ps.pitch_angle)


    def propagate_substep(self, propagator, substep):
        """Propagate substep in the propagator class and update the substep property of the particle state.
        """
        self.ps.substep = substep
        self.ps = propagator.move_substep(self.ps)
