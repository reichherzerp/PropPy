import pytest
import numpy as np
import pandas as pd
import os
#os.chdir('..')
#from proppy import proppy as pp
import proppy as pp


def test_point_source_oriented():
    print('\n----------------------------------')
    print('-> unit_test_point_source_oriented')
    energy = 10**10 #eV
    pos = [0,0,0]
    nr_particles = 10**2
    pitch_angle = 2*np.pi * 54.74/360
    phi = np.pi/4.0
    point_source_oriented = pp.PointSourceOriented(energy, pos, nr_particles, pitch_angle, phi)
    # test if the number of generated particles is correct
    assert nr_particles == len(point_source_oriented.particles)
    # test if the energy of a generated particle is correct
    assert energy == point_source_oriented.particles[0].ps.energy
    # check the position of the particles
    expected_pos = list(np.array(pos, dtype=np.float32))
    actual_pos = point_source_oriented.particles[0].ps.pos.tolist()
    assert expected_pos == actual_pos
    # check the direction of the particles
    precision = 4
    expected_direction = list(np.around(np.array([
        np.cos(phi)*np.sin(pitch_angle), 
        np.sin(phi)*np.sin(pitch_angle), 
        np.cos(pitch_angle)
    ], dtype=np.float32), precision))
    actual_direction = list(np.around(np.array(point_source_oriented.particles[0].ps.dir), precision))
    assert expected_direction == actual_direction