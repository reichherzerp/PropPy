from numba import float32, int32
import numpy as np
from numba.experimental import jitclass
from .magnetic_field import *
from .particle_state import *
from abc import ABCMeta, abstractmethod


constants_spec = [
    ('unit_distance', int32),
    ('speed', float32)
]

@jitclass(constants_spec)
class Constants():
    def __init__(self, unit_distance = 0):
        self.unit_distance = unit_distance
        if unit_distance == 0:
            self.speed = 2.998*10**8 # speed of light [m/s]
        elif unit_distance == 1:
            self.speed = 2.998*10**5 # speed of light [km/s]
        elif unit_distance == 2:
            self.speed = 9.716e-9 # speed of light [pc/s]

    def print_unit_kappa(self):
        if self.unit_distance == 0:
            return 'm²/s'
        elif self.unit_distance == 1:
            return 'km²/s'
        elif self.unit_distance == 2:
            return 'pc²/s'
        else:
            return ''

    def print_unit_speed(self):
        if self.unit_distance == 0:
            return 'm/s'
        elif self.unit_distance == 1:
            return 'km/s'
        elif self.unit_distance == 2:
            return 'pc/s'
        else:
            return ''

    def print_unit_distance(self):
        if self.unit_distance == 0:
            return 'm'
        elif self.unit_distance == 1:
            return 'km'
        elif self.unit_distance == 2:
            return 'pc'
        else:
            return ''