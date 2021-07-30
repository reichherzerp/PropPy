import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm

class Statistics():
    def __init__(self, df, dimensions):
        print('init statistics plotting class')
        self.df = df
        self.dimensions = dimensions
    
    def plot_distribution(self, axis, step, bins, file_name):
        # if user wants to plot the laast step, deduce the correct step number of the last step
        if step == -1:
            # get last step
            step = max(self.df['i'].tolist())
        # filter data according to the selected step
        df_i = self.df[self.df['i'] == step]
        
        # fit a normal distribution to the data:
        mu, std = norm.fit(df_i[axis])
        
        plt.figure(figsize=(4,4))
        # compute the factor between the original and normalized hists
        x_norm, b, p = plt.hist(df_i[axis], bins, density=True)
        x_original, b, p = plt.hist(df_i[axis], bins)
        scale = max(x_original)/max(x_norm)
    
        # plot the PDF scaled with the factor between the original and normalized hists
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)*scale
        plt.plot(x, p, 'k', linewidth=2)
        
        plt.xlabel(axis + ' [m]')
        plt.ylabel('# particles')
        label = "mu = %.2f,  std = %.2f" % (mu, std)
        #plt.legend()
        if file_name is not None:
            plt.tight_layout()
            plt.savefig(file_name)
        plt.show()
        
    def plot_diffusion_coefficients(self):
        nr_particles = len(list(map(int, (set(self.df['id'])))))
        df = self.df.sort_values('d')
        x = df['x'].values
        y = df['y'].values
        z = df['z'].values
        t = df['d'].values
        times = []
        kappa_xx = []
        kappa_yy = []
        kappa_perp = []
        kappa_zz = []
        for j in range(1, int(len(x)/nr_particles)):
            t_j = t[j*nr_particles]
            kappa_xx_current = 0
            kappa_yy_current = 0
            kappa_zz_current = 0
            for i in range(nr_particles):
                x_i = x[j*nr_particles+i]
                y_i = y[j*nr_particles+i]
                z_i = z[j*nr_particles+i]
                kappa_xx_current = kappa_xx_current + x_i**2
                kappa_yy_current = kappa_yy_current + y_i**2
                kappa_zz_current = kappa_zz_current + z_i**2
            kappa_xx.append(kappa_xx_current/(2*t_j))
            kappa_yy.append(kappa_yy_current/(2*t_j))
            kappa_perp.append((kappa_xx_current+kappa_yy_current)/(4*t_j))
            kappa_zz.append(kappa_zz_current/(2*t_j))
            times.append(t_j)
        plt.figure(figsize=(4,4))
        plt.plot(times, kappa_xx, label='$\kappa_{xx}$')
        plt.plot(times, kappa_yy, label='$\kappa_{yy}$')
        plt.plot(times, kappa_perp, label='$\kappa_\perp$')
        plt.plot(times, kappa_zz, label='$\kappa_\parallel$')
        plt.loglog()
        plt.legend()
        plt.show()
        
        d = {'kappa_perp':kappa_perp,'kappa_para':kappa_zz}
        return pd.DataFrame(d)