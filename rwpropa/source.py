"""The Source initializes the particles in the beginning of the simulation.

There are different sources available that can be customized by the user. 
The source specifies the initial state of the particles.

    Typical usage example:

    import rwpropa as rw
    
    nr_particles = 1*10**3
    source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    delta_rho_div_phi = 1 #1/2**0.5 # (delta_r_rho / delta_r_phi)
    energy = 3*10**15 # eV
    phi = 0.0
    pitch_angle = 2*np.pi * 54.74/360 # pitch angle for equal components in all directions

    sim = rw.Simulation()
    source = rw.PointSourceOriented(energy, source_pos, nr_particles, pitch_angle, phi)
    sim.add_source(source)
    sim.source.get_description()
"""


from .particle import Particle
import numpy as np
from abc import ABCMeta, abstractmethod


class SourceMeta(ABCMeta):
    """ Abstract meta class to check if all required attributes are implemented in the 
    sub classes.
    """

    required_attributes = []

    def __call__(self, *args, **kwargs):
        """ Checks if required attributes that have to be implemented in __init__ of all
        sub classes are really implemented. 

        Raises:
            ValueError: an error if not all required attributes are implemented.
        """
        obj = super(SourceMeta, self).__call__(*args, **kwargs)
        for attr_name in obj.required_attributes:
            if getattr(obj, attr_name) is None:
                raise ValueError('required attribute (%s) not set' % attr_name)
        return obj


class Source(object, metaclass=SourceMeta):
    """Abstract base class for all special sources.
    
    Functions with the label @abstractmethod have to be implemented in the special 
    observer classes.

    Attributes:
        energy: Initial energy of the particle.
        nr_particles: Number of particles emitted from the source.
        dimensions: Number of dimensions.
        particles: Particles with the initial particle state.
        source: Special source.
    """

    required_attributes = [
        'energy', 
        'nr_particles', 
        'dimensions', 
        'particles'
    ]  

    @abstractmethod
    def __init__(self, order):
        """Implementation required in all sub classes. All required_attributes have to 
        be implemented in sub classes.
        """
        pass


    def init_source(self):
        """Initialize parameters that are common for all special source classes.
        """
        self.particles = []
        self.dimensions = 3
        self.pos = np.array(self.pos, dtype=np.float32)
        self.inject_particles()

    
    def set_dimensions(self, dimensions):
        self.dimensions = dimensions


    @abstractmethod
    def inject_particles(self):
        """Each sub class (special source) has to implement how particles should be injected. Here,
        self.particles will be filled.
        """
        pass


    def reset_source(self):
        """Reset source after simulation in order to repeat the simulation afterwards.
        """
        self.empty_source()
        self.inject_particles()


    def empty_source(self):
        """Remove all particles in the source.
        """
        self.particles = []


    @abstractmethod
    def get_description_source_type(sefl):
        pass


    def get_description(self):
        """Print the information of the relevant parameters and the description of 
        the special source type that was chosen
        """
        self.get_description_general()
        self.get_description_parameters()
        self.get_description_source_type()


    def get_description_general(self):
        """Called by all special observer classes below.
        Introduction of the description output.
        """
        print("""Description Source:
                The source defines the start conditions of the particles 
                and covers the position, direction, energy, etc\n""")


    def get_description_parameters(self):   
        """Called by all special observer classes below.
        print out all relevant instance parameters.
        """
        print('position: ' , self.pos)
        print('number particles: ' ,self.nr_particles)
        print('energy: ' ,self.energy, ' eV')
         


class PointSourceOriented(Source):
    """A point source that emits particles into a user-defined direction.

    All particles start from a single point defined by the source position in 
    the user-defined direction. All particles have the exact same state in the 
    beginning.

    Attributes:
        energy: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        pos: A list that specify the source position. 
        nr_particles: An int that defines how many particles should be emitted.
        pitch_angle: The initial pitch angle of the particle -> angle between B-field and particle direction.
        phi: The initial angle between the particle direction in the xy-plane and the x-axis.
        particles: List of particles in the source. This list will be used in the simulation.
    """

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



class PointSourceIsotropicPhi(Source):
    """A point source that emits particles isotropically in phi.

    All particles start from a single point defined by the source position in 
    the user-defined direction. All particles have the exact same state in the 
    beginning, except for the direction, which is isotropic in phi.

    Attributes:
        energy: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        pos: A list that specify the source position. 
        nr_particles: An int that defines how many particles should be emited.
        particles: List of particles in the source. This list will be used in the simulation.
    """

    def __init__(self, energy, pos, nr_particles):
        self.energy = energy
        self.pos = pos
        self.nr_particles = nr_particles
        self.init_source()


    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            offset = 0
            if np.random.random() > 0.5:
                offset = 180
            pitch_angle = 2*np.pi * (54.74+offset)/360 # pitch angle for equal components in all directions
            phi = np.random.random()*2*np.pi
            particle_id = i
            p = Particle(particle_id, self.energy, self.pos[:], phi, pitch_angle, self.dimensions)
            self.particles.append(p)


    def get_description_source_type(self):
        print('source tpye: PointSourceIsotropic')



class PointSourceIsotropic(Source):
    """A point source that emits particles isotropically.

    All particles start from a single point defined by the source position in 
    the user-defined direction. All particles have the exact same state in the 
    beginning, except for the direction, which is isotropic.

    Attributes:
        energy: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        pos: A list that specify the source position. 
        nr_particles: An int that defines how many particles should be emited.
        particles: List of particles in the source. This list will be used in the simulation.
    """

    def __init__(self, energy, pos, nr_particles):
        self.energy = energy
        self.pos = pos
        self.nr_particles = nr_particles
        self.init_source()


    def sample_isotropic_vecotrs(self):
        """Samlpe correct isotropic vectors.

        See discussion in https://mathworld.wolfram.com/SpherePointPicking.html.
        """
        phi = np.random.random()*2*np.pi
        cos_pitch_angle = np.random.random()*2-1
        pitch_angle = np.arccos(cos_pitch_angle)
        return phi, pitch_angle


    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            phi, pitch_angle = self.sample_isotropic_vecotrs()
            particle_id = i
            p = Particle(particle_id, self.energy, self.pos[:], phi, pitch_angle, self.dimensions)
            self.particles.append(p)


    def get_description_source_type(self):
        print('source tpye: PointSourceIsotropic')



class SphereSourceIsotropic(Source):
    """A spherical source that emits particles isotropically.

    All particles start from a random point within a user defined sphere. 
    All particles have the exact same state in the 
    beginning, except for the position and the direction, which is isotropic.

    Attributes:
        energy: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        pos: A list that specify the source position. 
        nr_particles: An int that defines how many particles should be emited.
        particles: List of particles in the source. This list will be used in the simulation.
    """

    def __init__(self, energy, pos, nr_particles, radius):
        self.energy = energy
        self.pos = pos
        self.nr_particles = nr_particles
        self.radius = radius
        self.init_source()


    def sample_isotropic_vecotrs(self):
        """Samlpe correct isotropic vectors.

        See discussion in https://mathworld.wolfram.com/SpherePointPicking.html.
        """
        phi = np.random.random()*2*np.pi
        cos_pitch_angle = np.random.random()*2-1
        pitch_angle = np.arccos(cos_pitch_angle)
        rho = np.random.random()
        return phi, pitch_angle, rho


    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            phi, pitch_angle, rho = self.sample_isotropic_vecotrs()
            phi_pos, pitch_angle_pos, rho_pos = self.sample_isotropic_vecotrs()
            particle_id = i
            r = rho_pos**(1/3.0)
            pos = np.array([r*self.radius*np.cos(phi_pos)*np.sin(pitch_angle_pos), r*self.radius*np.sin(phi_pos)*np.sin(pitch_angle_pos), r*self.radius*np.cos(pitch_angle_pos)], dtype=np.float32)
            p = Particle(particle_id, self.energy, pos, phi, pitch_angle, self.dimensions)
            #p.set_random_direction(self, np.array([-1.0, 1.0, 1.0], dtype=np.float32))
            self.particles.append(p)


    def get_description_source_type(self):
        print('source tpye: PointSourceIsotropic')


class SphereSourceSurfaceIsotropic(Source):
    """A spherical source that emits particles isotropically from its surface.

    All particles start from a random point from the surface of a user defined sphere. 
    All particles have the exact same state in the 
    beginning, except for the position and the direction, which is isotropic.

    Attributes:
        energy: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        pos: A list that specify the source position. 
        nr_particles: An int that defines how many particles should be emited.
        particles: List of particles in the source. This list will be used in the simulation.
    """

    def __init__(self, energy, pos, nr_particles, radius):
        self.energy = energy
        self.pos = pos
        self.nr_particles = nr_particles
        self.radius = radius
        self.init_source()


    def sample_isotropic_vecotrs(self):
        """Samlpe correct isotropic vectors.

        See discussion in https://mathworld.wolfram.com/SpherePointPicking.html.
        """
        phi = np.random.random()*2*np.pi
        cos_pitch_angle = np.random.random()*2-1
        pitch_angle = np.arccos(cos_pitch_angle)
        rho = np.random.random()
        return phi, pitch_angle, rho


    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            phi, pitch_angle, rho = self.sample_isotropic_vecotrs()
            phi_pos, pitch_angle_pos, rho_pos = self.sample_isotropic_vecotrs()
            particle_id = i
            pos = np.array([self.radius*np.cos(phi_pos)*np.sin(pitch_angle_pos), self.radius*np.sin(phi_pos)*np.sin(pitch_angle_pos), self.radius*np.cos(pitch_angle_pos)], dtype=np.float32)
            p = Particle(particle_id, self.energy, pos, phi, pitch_angle, self.dimensions)
            self.particles.append(p)


    def get_description_source_type(self):
        print('source tpye: PointSourceIsotropic')



class QubicSourceIsotropic(Source):
    """A quibic source that emits particles isotropically from within its volume.

    All particles start from a random point from the 3d cube of a user defined width. 
    All particles have the exact same state in the 
    beginning, except for the position and the direction, which is isotropic.

    Attributes:
        energy: An b array specifying observed substeps [1_substep,2_substep,3_substep].
                  Only observing once per step: substeps = [False, False, True].
        pos: A list that specify the source position. 
        nr_particles: An int that defines how many particles should be emited.
        particles: List of particles in the source. This list will be used in the simulation.
    """

    def __init__(self, energy, pos, nr_particles, width):
        self.energy = energy
        self.pos = pos
        self.nr_particles = nr_particles
        self.width = width
        self.init_source()

    def sample_isotropic_vecotrs(self):
        """Samlpe correct isotropic vectors.

        See discussion in https://mathworld.wolfram.com/SpherePointPicking.html.
        """
        phi = np.random.random()*2*np.pi
        cos_pitch_angle = np.random.random()*2-1
        pitch_angle = np.arccos(cos_pitch_angle)
        rho = np.random.random()
        return phi, pitch_angle, rho

    def inject_particles(self):
        self.particles = []
        for i in range(self.nr_particles):
            particle_id = i
            phi, pitch_angle, rho = self.sample_isotropic_vecotrs()
            pos = np.array([(np.random.random()*2-1)*self.width, (np.random.random()*2-1)*self.width, (np.random.random()*2-1)*self.width], dtype=np.float32)
            p = Particle(particle_id, self.energy, pos, phi, pitch_angle, self.dimensions)
            self.particles.append(p)


    def get_description_source_type(self):
        print('source tpye: PointSourceIsotropic')