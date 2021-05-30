import matplotlib.pyplot as plt
import random
import numpy as np
from numba import jit, int32, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass

class Source():
    def __init__(self, nr_particles, position):
        self.nr_particles = nr_particles
        self.position = position