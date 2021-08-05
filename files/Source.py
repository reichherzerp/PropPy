from files.Particle import Particle
import numpy as np
from abc import ABC, abstractmethod


class Source(ABC):
    def __init__(self):
        super().__init__()
        print('Source initialized')

    def init_source(self):
        self.particles = []
        self.dimensions = 3

    @abstractmethod
    def inject_particles(self):
        pass

    def empty_source(self):
        self.particles = []
    


class PointSource(Source):
    def __init__(self, gyro_radius, pos, nr_particles):
        self.gyro_radius = gyro_radius
        self.pos = pos
        self.nr_particles = nr_particles
        self.init_source()

    def inject_particles(self):
        self.particles = []
        phi = 0.0
        pitch_angle = 2*np.pi * 54.74/360 # pitch angle for equal components in all directions
        for i in range(self.nr_particles):
            particle_id = i
            p = Particle(particle_id, self.gyro_radius, self.pos[:], phi, pitch_angle, self.dimensions)
            self.particles.append(p)