import unittest
from .observer import *
from .propagator import *
from .source import *

class TestSource(unittest.TestCase):

    def test_point_source_oriented(self):
        energy = 10**10 #eV
        pos = [0,0,0]
        nr_particles = 10**2
        pitch_angle = 1
        phi = 0
        point_source_oriented = PointSourceOriented(energy, pos, nr_particles, pitch_angle, phi)