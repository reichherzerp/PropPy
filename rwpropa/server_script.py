import time
import numpy as np
import pandas as pd
import os
os.chdir('..')
import rwpropa as rw
import sys


def simulate(file_name):
    sim = rw.PlasmoidSimulation()
    start = time.time()
    sim.simulate(file_name)
    print('finished: ', file_name)
    end = time.time()
    print("Elapsed (after compilation) = %s" % (end - start))


if __name__ == '__main__':
    simulate('data_'+str(sys.argv[1]))