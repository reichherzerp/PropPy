import time
import numpy as np
import pandas as pd
import os
os.chdir('..')
import rwpropa as rw
import sys


def simulate(file_name, i):
    file_name_i = file_name + str(i)
    parameters = {'nr_particles': 10**2, 'radius': 10**14 }
    nr_particles = 10**2 
    radius = 10**14 
    energy = 10**15
    nr_steps = 1*10**4
    diffusion_coefficient = 1.5*10**20
    step_size = 1*10**12
    sim = rw.PlasmoidSimulation(nr_particles = nr_particles, radius = radius, energy = energy, nr_steps = nr_steps, diffusion_coefficient = diffusion_coefficient, step_size = step_size)
    start = time.time()
    sim.simulate(file_name_i)
    print('finished: ', i)
    end = time.time()
    print('Elapsed (after compilation) = %s' % (end - start))

    # open file for writing
    f = open(file_name+'_parameters.txt','w')
    # write file
    f.write( str(parameters) )
    # close file
    f.close()


if __name__ == '__main__':
    file_name = str(sys.argv[1])
    i = str(sys.argv[2])
    simulate(file_name, i)