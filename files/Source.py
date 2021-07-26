from files.Particle import Particle

class Source():
    def __init__(self, gyro_radius, mean_free_path, pos, nr_particles):
        print('source')
        self.particles = []
        for i in range(nr_particles):
            particle_id = i
            p = Particle(particle_id, gyro_radius, mean_free_path, pos)
            self.particles.append(p)