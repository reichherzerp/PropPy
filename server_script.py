import time
import proppy.simulation as simulation
import sys
import os
os.chdir('..')
import proppy as pp


def simulate(file_name, i, total_number_simulations):
    file_name_i = file_name + str(i)
    param = {
        'nr_particles': 10**2, 
        'radius': 10**14,
        'energy': 10**15,
        'nr_steps': 1*10**4,
        'diffusion_coefficient': 1.5*10**20,
        'step_size': 1*10**12,
        'total_number_simulations': total_number_simulations
    }
    
    sim = rw.PlasmoidSimulation(nr_particles = param['nr_particles'], radius = param['radius'], energy = param['energy'], nr_steps = param['nr_steps'], diffusion_coefficient = param['diffusion_coefficient'], step_size = param['step_size'])
    start = time.time()
    sim.simulate(file_name_i)
    print('finished: ', i)
    end = time.time()
    print('Elapsed (after compilation) = %s' % (end - start))

    if i == 1: #only save parameters once
        # open file for writing
        f = open(file_name+'_parameters.txt','w')
        # write file
        f.write( str(param) )
        # close file
        f.close()


if __name__ == '__main__':
    # get parameters from batch script (job.sh)
    file_name = str(sys.argv[1])
    # current sub simulation id
    i = str(sys.argv[2])
    # number of sub simulations
    total_number_simulations = str(sys.argv[3])
    # run sub simulation 
    simulate(file_name, i, total_number_simulations)