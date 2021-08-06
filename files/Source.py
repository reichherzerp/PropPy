from files.Particle import Particle
import numpy as np
from abc import ABCMeta, abstractmethod


class SourceMeta(ABCMeta):
    # required attributes that have to be implemented in __init__ of all
    # sub classes
    required_attributes = []

    def __call__(self, *args, **kwargs):
        # check if required attributes that have to be implemented in __init__ of all
        # sub classes are really implemented. Raise an error if not
        obj = super(SourceMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj


class Source(object, metaclass=SourceMeta):
    # all required_attributes have to be implemented in sub classes
    required_attributes = ['energy', 'nr_particles', 'dimensions', 'particles']  

    @abstractmethod
    def __init__(self, order):
        # implementation required in all sub classes.
        # all required_attributes have to be implemented in sub classes
        pass

    def init_source(self):
        # initialize parameters that are common for all special source classes
        self.particles = []
        self.dimensions = 3
        self.inject_particles()

    @abstractmethod
    def inject_particles(self):
        # each sub class (special source) has to implement how particles should be injected. Here,
        # self.particles will be filled
        pass

    def reset_source(self):
        # reset source after simulation in order to repeat the simulation afterwards
        self.empty_source()
        self.inject_particles()

    def empty_source(self):
        # remove all particles in the source
        self.particles = []
    


class PointSourceOriented(Source):
    def __init__(self, gyro_radius, pos, nr_particles, pitch_angle, phi):
        self.energy = 10
        self.gyro_radius = gyro_radius
        self.pos = pos
        self.nr_particles = nr_particles
        self.pitch_angle = pitch_angle
        self.phi = phi
        self.init_source()

    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            particle_id = i
            p = Particle(particle_id, self.gyro_radius, self.pos[:], self.phi, self.pitch_angle, self.dimensions)
            self.particles.append(p)



class PointSourceIsotropic(Source):
    def __init__(self, gyro_radius, pos, nr_particles):
        self.energy = 10
        self.gyro_radius = gyro_radius
        self.pos = pos
        self.nr_particles = nr_particles
        self.init_source()

    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            phi = np.random * 360
            pitch_angle = np.random * 2*np.pi
            particle_id = i
            p = Particle(particle_id, self.gyro_radius, self.pos[:], phi, pitch_angle, self.dimensions)
            self.particles.append(p)