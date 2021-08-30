import unittest
import os
os.chdir('..')
import rwpropa as rw

class TestSource(unittest.TestCase):

    def test_point_source_oriented(self):
        print('-> test_point_source_oriented')
        energy = 10**10 #eV
        pos = [0,0,0]
        nr_particles = 10**2
        pitch_angle = 0
        phi = 0
        point_source_oriented = rw.PointSourceOriented(energy, pos, nr_particles, pitch_angle, phi)
        # test if the number of generated particles is correct
        self.assertEqual(nr_particles, len(point_source_oriented.particles))
        # test if the energy of a generated particle is correct
        self.assertEqual(energy, point_source_oriented.particles[0].ps.energy)
        

    