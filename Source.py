import matplotlib.pyplot as plt
import random
import numpy as np
from numba import jit, int32, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass

class Source():
    # generate particles
    @jit(nopython=True)
    def source():
        particles = List()
        for j in range(nr_particles):
            direction = List()
            random_values = np.random.randint(low=0, high=1, size=3)*2-1
            [direction.append(x) for x in random_values]
            particles.append(particle(gyro_radius, free_mean_path_para, free_mean_path_perp, direction))
        return particles