import numpy as np
import pandas as pd
import os
os.chdir('../..')
import rwpropa as rw

sim = rw.Simulation()

nr_particles = 10**3
source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
energy = 3*10**15 # eV

source = rw.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
sim.add_source(source)
sim.source.get_description()