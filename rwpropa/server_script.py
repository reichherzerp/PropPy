import time
import numpy as np
import pandas as pd
import os
os.chdir('..')
import rwpropa as rw
import sys


def simulate(file_name):
    nr_particles = 10**2 
    radius = 10**14 
    energy = 10**15
    nr_steps = 1*10**4
    diffusion_coefficient = 1.5*10**20
    step_size = 1*10**12
    sim = rw.PlasmoidSimulation(nr_particles = nr_particles, radius = radius, energy = energy, nr_steps = nr_steps, diffusion_coefficient = diffusion_coefficient, step_size = step_size)
    start = time.time()
    sim.simulate(file_name)
    print('finished: ', file_name)
    end = time.time()
    print("Elapsed (after compilation) = %s" % (end - start))


if __name__ == '__main__':
    file_name = str(sys.argv[1])
    simulate(file_name)