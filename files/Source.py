from files.Particle import Particle
import numpy as np
from abc import ABC, ABCMeta, abstractmethod


class SourceMeta(ABCMeta):
    # required attributes that have to be implemented in __init__ of all
    # sub classes
    required_attributes = []

    def __call__(self, *args, **kwargs):
        # check if required attributes that have to be implemented in __init__ of all
        # sub classes are really implemented. Raise an error if not
        obj = super(SourceMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if not getattr(obj, attr_name):
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj


class Source(object, metaclass=SourceMeta):
    # all required_attributes have to be implemented in sub classes
    required_attributes = ['energy', 'nr_particles']

    @abstractmethod
    def __init__(self, order):
        # implementation required in all sub classes.
        # all required_attributes have to be implemented in sub classes
        pass

    def init_source(self):
        # initialize parameters that are common for all special source classes
        self.particles = []
        self.dimensions = 3

    @abstractmethod
    def inject_particles(self):
        #
        pass

    def empty_source(self):
        # remove all particles in the source
        self.particles = []
    


class PointSource(Source):
    def __init__(self, gyro_radius, pos, nr_particles):
        self.energy = 10
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