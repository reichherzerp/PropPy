import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
os.chdir('..')
import rwpropa as rw

trajectory_lengths = np.load('merged.npy')
diff = 1.5*10**20
radius = 10**14
nr_particles = 10**6
plt.figure(figsize=(5,3))
bins = 30
d = trajectory_lengths/10**14
hist, bins = np.histogram(d, bins=bins)
logbins = np.logspace(np.log10(min(d)),np.log10(max(d)),len(bins))
plt.hist(d, bins=logbins, alpha=0.5, label='$\kappa =$ {:.1e}m$^2$/s'.format(diff))
plt.axvline(x=1, color='k', ls='--', label='plasmoid radius')
plt.title('total # particles = {:.0e}'.format(len(d)))
plt.xlabel('D/{:.0e}m'.format(radius))
plt.ylabel('# particles')
plt.loglog()
plt.legend()
plt.show()