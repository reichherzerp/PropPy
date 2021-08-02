from files.Particle import Particle

class Source():
    def __init__(self, gyro_radius, mean_free_path, pos, nr_particles):
        self.gyro_radius = gyro_radius
        self.mean_free_path = mean_free_path
        self.pos = pos
        self.nr_particles = nr_particles
        print('source')
        self.particles = []

    def init_source(self):
        self.particles = []
        for i in range(self.nr_particles):
            particle_id = i
            p = Particle(particle_id, self.gyro_radius, self.mean_free_path, self.pos[:])
            self.particles.append(p)