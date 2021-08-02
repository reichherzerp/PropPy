from files.Particle import Particle

class Source():
    def __init__(self, gyro_radius, pos, nr_particles, dimensions):
        self.gyro_radius = gyro_radius
        self.pos = pos
        self.nr_particles = nr_particles
        self.dimensions = dimensions
        self.particles = []
        print('Source initialized')

    def init_source(self):
        self.particles = []
        for i in range(self.nr_particles):
            particle_id = i
            p = Particle(particle_id, self.gyro_radius, self.pos[:], self.dimensions)
            self.particles.append(p)