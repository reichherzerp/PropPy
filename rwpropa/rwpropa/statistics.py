import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm

class Statistics():    
    def __init__(self, df):
        print('init statistics plotting class')
        self.df = df
        self.dimensions = 3


    def set_dimensions(self, dimensions):
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
        
        xlabel = axis
        if axis == 'x' or axis == 'y' or axis == 'z' or axis == 'd':
            xlabel = xlabel + ' [m]'
        plt.xlabel(xlabel)
        plt.ylabel('# particles')
        label = "mu = %.2f,  std = %.2f" % (mu, std)
        #plt.legend()
        if file_name is not None:
            plt.tight_layout()
            plt.savefig(file_name)
        plt.show()

        
    def plot_diffusion_coefficients(self, error):
        nr_particles = len(list(map(int, (set(self.df['id'])))))
        # sort the pandas dataframe so that the df can be distributed in packages of length=nr_particles
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
        if error:
            kappa_xx_err = []
            kappa_yy_err = []
            kappa_perp_err = []
            kappa_zz_err = []
        c = 299792458 # speed of light [m/s]
        nr_steps =  int(len(d)/nr_particles)
        # calculate the running diffusion coefficient kappa_i(t) for each step
        # running diffusion coefficients: kappa_i(t) = <x_i>^2/(2t) = <x_i>^2/(2d/c)
        for j in range(nr_steps):
            d_j = d[j*nr_particles]
            # get list of x_i^2 for all particles at the current step
            x_squared = np.array(x[j*nr_particles : (1+j)*nr_particles])**2
            y_squared = np.array(y[j*nr_particles : (1+j)*nr_particles])**2
            z_squared = np.array(z[j*nr_particles : (1+j)*nr_particles])**2
            # compute the running diffusion coefficient at t: kappa_i(t) = <x_i>^2/(2d/c)
            kappa_xx_running = np.mean(x_squared)/(2*d_j/c)
            kappa_yy_running = np.mean(y_squared)/(2*d_j/c)
            kappa_zz_running = np.mean(z_squared)/(2*d_j/c)
            kappa_xx.append(kappa_xx_running)
            kappa_yy.append(kappa_yy_running)
            kappa_perp.append((kappa_xx_running+kappa_yy_running)/2)
            kappa_zz.append(kappa_zz_running)
            times.append(d_j)
            if error:
                # compute the standard deviation of the mean square displacement in each step
                kappa_xx_running_err = np.std(x_squared)/(2*d_j/c)
                kappa_yy_running_err = np.std(y_squared)/(2*d_j/c)
                kappa_zz_running_err = np.std(z_squared)/(2*d_j/c)
                kappa_xx_err.append(kappa_xx_running_err)
                kappa_yy_err.append(kappa_yy_running_err)
                kappa_perp_err.append((kappa_xx_running_err+kappa_yy_running_err)/2)
                kappa_zz_err.append(kappa_zz_running_err)
        plt.figure(figsize=(4,4))
        if error:
            # plot with error bars
            plt.errorbar(times, kappa_zz, yerr=kappa_zz_err, fmt=".", elinewidth=0.5, markersize=4, c='dodgerblue', label='$\kappa_\parallel$')
            plt.errorbar(times, kappa_perp, yerr=kappa_perp_err, fmt=".", elinewidth=0.5, markersize=4, c='brown', label='$\kappa_\perp$')
        else:
            # plot without error bars
            plt.plot(times, kappa_zz, zorder=1, label='$\kappa_\parallel$', c='dodgerblue') 
            plt.scatter(times, kappa_zz, zorder=2, s=2, c='dodgerblue')
            plt.plot(times, kappa_perp, zorder=1, label='$\kappa_\perp$', c='brown') 
            plt.scatter(times, kappa_perp, zorder=2, s=2, c='brown')
        
        plt.xlabel('d [m]')
        plt.ylabel('running diffusion coefficients')
        plt.loglog()
        plt.legend()
        plt.show()
        
        d = {'d': times, 'kappa_perp':kappa_perp,'kappa_para':kappa_zz}
        return pd.DataFrame(d)

