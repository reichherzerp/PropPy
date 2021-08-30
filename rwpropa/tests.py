import unittest
import os
os.chdir('..')
import rwpropa as rw

class TestSource(unittest.TestCase):

    def test_point_source_oriented(self):
        energy = 10**10 #eV
        pos = [0,0,0]
        nr_particles = 10**2
        pitch_angle = 1
        phi = 0
        point_source_oriented = rw.PointSourceOriented(energy, pos, nr_particles, pitch_angle, phi)
        #self.assertEqual(
        self.assertEqual(energy, point_source_oriented.particles[0].ps.energy)
        #len(point_source_oriented.particles)