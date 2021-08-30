from .particle import Particle
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
    required_attributes = [
        'energy', 
        'nr_particles', 
        'dimensions', 
        'particles'
    ]  

    @abstractmethod
    def __init__(self, order):
        # implementation required in all sub classes.
        # all required_attributes have to be implemented in sub classes
        pass


    def init_source(self):
        # initialize parameters that are common for all special source classes
        self.particles = []
        self.dimensions = 3
        self.pos = np.array(self.pos, dtype=np.float32)
        self.inject_particles()

    
    def set_dimensions(self, dimensions):
        self.dimensions = dimensions


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


    @abstractmethod
    def get_description_source_type(sefl):
        pass


    def get_description(self):
        # print the information of the relevant parameters and the description of 
        # the special source type that was chosen
        self.get_description_general()
        self.get_description_parameters()
        self.get_description_source_type()


    def get_description_general(self):
        # called by all special observer classes below.
        # introduction of the description output
        print("""Description Source:
                The source defines the start conditions of the particles 
                and covers the position, direction, energy, etc\n""")


    def get_description_parameters(self):   
        # called by all special observer classes below.
        # print out all relevant instance parameters
        print('position: ' , self.pos)
        print('number particles: ' ,self.nr_particles)
        print('energy: ' ,self.energy, ' eV')
         


class PointSourceOriented(Source):
    def __init__(self, energy, pos, nr_particles, pitch_angle, phi):
        self.energy = energy
        self.pos = pos
        self.nr_particles = nr_particles
        self.pitch_angle = pitch_angle
        self.phi = phi
        self.init_source()


    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            particle_id = i
            p = Particle(particle_id, self.energy, self.pos[:], self.phi, self.pitch_angle, self.dimensions)
            self.particles.append(p)


    def get_description_source_type(self):
        print('source tpye: PointSourceOriented')
        print('pitch angle: ' ,self.pitch_angle)
        print('phi: ' ,self.phi)



class PointSourceIsotropic(Source):
    def __init__(self, energy, pos, nr_particles):
        self.energy = energy
        self.pos = pos
        self.nr_particles = nr_particles
        self.init_source()


    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            phi = np.random.rand() * 360
            pitch_angle = np.random.rand() * 2*np.pi
            particle_id = i
            p = Particle(particle_id, self.energy, self.pos[:], phi, pitch_angle, self.dimensions)
            self.particles.append(p)


    def get_description_source_type(self):
        print('source tpye: PointSourceIsotropic')