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
        d = df['d'].values
        times = []
        kappa_xx = []
        kappa_yy = []
        kappa_perp = []
        kappa_zz = []
        c = 299792458 # speed of light [m/s]
        nr_steps =  int(len(d)/nr_particles)
        # calculate the running diffusion coefficient kappa_i(t) for each step
        # running diffusion coefficients: kappa_i(t) = <x_i>^2/(2t) = <x_i>^2/(2x_i/c)
        for j in range(1, nr_steps):
            d_j = d[j*nr_particles]
            x_squared = np.array(x[j*nr_particles : (1+j)*nr_particles])**2
            y_squared = np.array(y[j*nr_particles : (1+j)*nr_particles])**2
            z_squared = np.array(z[j*nr_particles : (1+j)*nr_particles])**2
            kappa_xx_running = np.mean(x_squared)/(2*d_j/c)
            kappa_yy_running = np.mean(y_squared)/(2*d_j/c)
            kappa_zz_running = np.mean(z_squared)/(2*d_j/c)
            kappa_xx_running_err = np.std(x_squared)/(2*d_j/c)
            kappa_yy_running_err = np.std(y_squared)/(2*d_j/c)
            kappa_zz_running_err = np.std(z_squared)/(2*d_j/c)
            x_squared = 0
            y_squared = 0
            z_squared = 0
            kappa_xx.append(kappa_xx_running)
            kappa_yy.append(kappa_yy_running)
            kappa_perp.append((kappa_xx_running+kappa_yy_running)/2)
            kappa_zz.append(kappa_zz_running)
            times.append(d_j)
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