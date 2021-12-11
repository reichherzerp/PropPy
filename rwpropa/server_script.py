import time
import numpy as np
import pandas as pd
import os
os.chdir('..')
import rwpropa as rw
import sys


def simulate(file_name):
    sim = rw.PlasmoidSimulation(nr_particles = 10**6, radius = 10**14, energy = 10**15, nr_steps = 1*10**6, diffusion_coefficient = 1.5*10**20, step_size = 1*10**11)
    start = time.time()
    sim.simulate(file_name)
    print('finished: ', file_name)
    end = time.time()
    print("Elapsed (after compilation) = %s" % (end - start))


if __name__ == '__main__':
    file_name = str(sys.argv[1])
    simulate(file_name)