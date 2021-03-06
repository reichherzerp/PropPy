"""Scripts to compute statistical quantities of the particles.

As the particles propagate via a random walk, statistical properties 
of many particles are interesting, such as the diffusion coefficients and
particle distributions. These quantities can be compared to analytical 
predictions.


    Typical usage example:

    import proppy as pp

    df = pd.read_pickle("data/data_sim.pkl")
    df_time_evolution_observer = df.loc[df['radius'] == -1.0]
    sta = pp.Statistics(df_time_evolution_observer)
    errors = False
    df_kappas = sta.plot_diffusion_coefficients(errors)
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm


class Statistics():   
    """Statistics class for visualizing statistical properties of particle transport.

    After loading the simulation data, distributions of particles and running diffusion
    coefficients can be plotted.

    Attributes:
        dimensions: An int for defining the dimensions.
        df: A pandas dataftrame with the simulation data.
    """

    def __init__(self, df):
        print('init statistics plotting class')
        self.df = df
        self.dimensions = 3


    def set_dimensions(self, dimensions):
        self.dimensions = dimensions
    

    def plot_distribution(self, column, step, bins, file_name):
        """Plotting particle distributions at a given step number.
        
        The column of the df that should be plotted can be choosen freely. 
        Possible is to plot the particle distribution along directions, but
        also all other dataframe (df) columns are possible and sometimes usefull.

        Args:
            column: An int that specifies which column to plot.
            step: An int that defines the step number -> row of the df.
            bins: An int that defines the number of bins in the plotted histogram.
            file_name: String or None. If not None, the plot will be saved with the given String as a name.
        """

        # if user wants to plot the laast step, deduce the correct step number of the last step
        if step == -1:
            # get last step
            step = max(self.df['i'].tolist())
        # filter data according to the selected step
        df_i = self.df[self.df['i'] == step]
        
        # fit a normal distribution to the data:
        mu, std = norm.fit(df_i[column])
        
        plt.figure(figsize=(4,4))
        # compute the factor between the original and normalized hists
        x_norm, b, p = plt.hist(df_i[column], bins, density=True)
        x_original, b, p = plt.hist(df_i[column], bins)
        scale = max(x_original)/max(x_norm)
    
        # plot the PDF scaled with the factor between the original and normalized hists
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)*scale
        plt.plot(x, p, 'k', linewidth=2)
        
        xlabel = column
        if column == 'x' or column == 'y' or column == 'z' or column == 'd':
            xlabel = xlabel + ' [m]'
        plt.xlabel(xlabel)
        plt.ylabel('# particles')
        label = "mu = %.2f,  std = %.2f" % (mu, std)
        #plt.legend()
        if file_name is not None:
            plt.tight_layout()
            plt.savefig(file_name)
        plt.show()


    def plot_arrival_times(self, radius, nr_particles, diffusion_coefficient, file_name, bins_log=20, bins_lin=50, d_max=None):
        """Plotting the histograms for arrival times on a sphere.
        
        Two plots are created that show the number of particles arriving the
        sphere as a function of time:
        - in lin-lin representation to resolve later times
        - in log-log representation to resolve initial times

        Args:
            radius: A double that defines the radius of the sphere.
            nr_particles: An int only relevant for the title of the plots.
            diffusion_coefficient: A float only relevant for the title of the plots.
            file_name: A string that points to the source of the data that should be plotted.
            bins_log: An int that defines the number of bins for the lin-lin hist.
            bins_lin: An int that defines the number of bins for the log-log hist.
            d_max: An float that defines the maximum trajectory length that is shown.
        """
        df = pd.read_pickle(file_name+'.pkl')
        trajectory_lengths = df['d']
        plt.figure(figsize=(5,3))
        bins = bins_log
        d = trajectory_lengths/10**14
        if d_max == None:
            d_max_value = max(d)
        else:
            d_max_value = d_max
        print(max(d))
        hist, bins = np.histogram(d, bins=bins)
        logbins = np.logspace(np.log10(min(d)),np.log10(d_max_value),len(bins))
        plt.hist(d, bins=logbins, alpha=0.5, label=r'$\kappa =$ {:.1e}m$^2$/s'.format(diffusion_coefficient))
        plt.axvline(x=1, color='k', ls='--', label='plasmoid radius')
        plt.title('total # particles = {:.0e}'.format(nr_particles))
        plt.xlabel('D/{:.0e}m'.format(radius))
        plt.ylabel('# particles')
        plt.loglog()
        plt.legend()
        plt.show()
        
        plt.figure(figsize=(5,3))
        bins = bins_lin
        d = trajectory_lengths/10**14
        hist, bins = np.histogram(d, bins=bins) 
        linbins = np.linspace((min(d)),(d_max_value),len(bins))
        plt.hist(d, bins=linbins, alpha=0.5, label=r'$\kappa =$ {:.1e}m$^2$/s'.format(diffusion_coefficient))
        plt.axvline(x=1, color='k', ls='--', label='plasmoid radius')
        plt.title('total # particles = {:.0e}'.format(nr_particles))
        plt.xlabel('D/{:.0e}m'.format(radius))
        plt.ylabel('# particles')
        plt.legend()
        plt.show()


    def get_diffusion_coefficients(self, speed_light=3*10**8):
        nr_particles = len(list(map(int, (set(self.df['id'])))))
        # sort the pandas dataframe so that the df can be distributed in packages of length=nr_particles
        df = self.df.sort_values('d')
        x = df['x'].values
        y = df['y'].values
        z = df['z'].values
        d = df['d'].values
        l = []
        kappa_xx = []
        kappa_yy = []
        kappa_perp = []
        kappa_zz = []
        kappa_xx_err = []
        kappa_yy_err = []
        kappa_perp_err = []
        kappa_zz_err = []
        kappa = []
        c = speed_light # speed of light
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
            kappa.append((kappa_xx_running+kappa_yy_running+kappa_zz_running)/3)
            l.append(d_j)
            # compute the standard deviation of the mean square displacement in each step
            kappa_xx_running_err = np.std(x_squared)/(2*d_j/c)
            kappa_yy_running_err = np.std(y_squared)/(2*d_j/c)
            kappa_zz_running_err = np.std(z_squared)/(2*d_j/c)
            kappa_xx_err.append(kappa_xx_running_err)
            kappa_yy_err.append(kappa_yy_running_err)
            kappa_perp_err.append((kappa_xx_running_err+kappa_yy_running_err)/2)
            kappa_zz_err.append(kappa_zz_running_err)
        d = {
            'l': l, 
            'kappa': kappa,
            'kappa_xx': kappa_xx,
            'kappa_xx_err': kappa_xx_err,
            'kappa_yy': kappa_yy,
            'kappa_yy_err': kappa_yy_err,
            'kappa_zz': kappa_zz,
            'kappa_zz_err': kappa_zz_err,
            'kappa_perp': kappa_perp,
            'kappa_perp_err': kappa_perp_err,
            'kappa_para': kappa_zz
        }
        return pd.DataFrame(d)

        
    def plot_diffusion_coefficients(self, isotropic=True, error=False, file_name=None, plot=True, n_points_plateau = 100):
        """Plotting the running diffusion coefficients.
        
        The computation is described in:
        [1]Reichherzer, P., Becker Tjus, J., Zweibel, E. G., Merten, L., and Pueschel, M. J., 
        ???Turbulence-level dependence of cosmic ray parallel diffusion???, 
        Monthly Notices of the Royal Astronomical Society, vol. 498, no. 4, pp. 5051???5064, 2020. 
        doi:10.1093/mnras/staa2533.

        Args:
            isotropic: A bool to define if the diffusion is isotropic.
            error: A bool to define if the error bars should be plotted.
        """

        df_kappas =  self.get_diffusion_coefficients()
        l = df_kappas['l'].values.tolist()
        kappa_xx = df_kappas['kappa_xx'].values.tolist()
        kappa_yy = df_kappas['kappa_yy'].values.tolist()
        kappa_zz = df_kappas['kappa_zz'].values.tolist()
        kappa_xx_err = df_kappas['kappa_xx_err'].values.tolist()
        kappa_yy_err = df_kappas['kappa_yy_err'].values.tolist()
        kappa_zz_err = df_kappas['kappa_zz_err'].values.tolist()
        kappa_perp = df_kappas['kappa_perp'].values.tolist()
        kappa_perp_err = df_kappas['kappa_perp_err'].values.tolist()
        plt.figure(figsize=(3,2))
        if error:
            # plot with error bars
            if isotropic:
                plt.errorbar(l, kappa_xx, yerr=kappa_xx_err, fmt=".", elinewidth=0.5, markersize=4, c='k', label=r'$\kappa_{xx}$')
                plt.errorbar(l, kappa_yy, yerr=kappa_yy_err, fmt=".", elinewidth=0.5, markersize=4, c='dodgerblue', label=r'$\kappa_{yy}$')
                plt.errorbar(l, kappa_zz, yerr=kappa_zz_err, fmt=".", elinewidth=0.5, markersize=4, c='brown', label=r'$\kappa_{zz}$')
            else:
                plt.errorbar(l, kappa_zz, yerr=kappa_zz_err, fmt=".", elinewidth=0.5, markersize=4, c='dodgerblue', label=r'$\kappa_\parallel$')
                plt.errorbar(l, kappa_perp, yerr=kappa_perp_err, fmt=".", elinewidth=0.5, markersize=4, c='brown', label=r'$\kappa_\perp$')
        else:
            # plot without error bars
            if isotropic:
                plt.plot(l, kappa_xx, zorder=1, label=r'$\kappa_{xx}$', c='k') 
                plt.scatter(l, kappa_xx, zorder=2, s=2, c='k')
                plt.plot(l, kappa_yy, zorder=1, label=r'$\kappa_{yy}$', c='brown') 
                plt.scatter(l, kappa_yy, zorder=2, s=2, c='brown')
                plt.plot(l, kappa_zz, zorder=1, label=r'$\kappa_{zz}$', c='dodgerblue') 
                plt.scatter(l, kappa_zz, zorder=2, s=2, c='dodgerblue')
            else:
                plt.plot(l, kappa_perp, zorder=1, label=r'$\kappa_\perp$', c='brown') 
                plt.scatter(l, kappa_perp, zorder=2, s=2, c='brown')
                plt.plot(l, kappa_zz, zorder=1, label=r'$\kappa_\parallel$', c='dodgerblue') 
                plt.scatter(l, kappa_zz, zorder=2, s=2, c='dodgerblue')
        
        plt.xlabel('trajectory length [m]')
        plt.ylabel('running diff. coeff. [m$^2$/s]')
        plt.loglog()
        plt.legend()
        if file_name is not None:
            plt.tight_layout()
            plt.savefig(file_name)
        if plot:
            plt.show()

        n = n_points_plateau
        print('diffusion coefficients computed between ' + str("{:.2e}".format(l[-n_points_plateau])) + 'm and ' + str("{:.2e}".format(l[-1])) +'m with ' + str(n) + ' data points')
        print('kappa_{xx}:', f"{np.mean(kappa_xx[-n:]):.3}", 'm??/s', '+-', f"{np.std(kappa_xx[-n:]):.3}", 'm??/s')
        print('kappa_{yy}:', f"{np.mean(kappa_yy[-n:]):.3}", 'm??/s', '+-', f"{np.std(kappa_yy[-n:]):.3}", 'm??/s')
        print('kappa_{zz}:', f"{np.mean(kappa_zz[-n:]):.3}", 'm??/s', '+-', f"{np.std(kappa_zz[-n:]):.3}", 'm??/s')
        
        d = {
            'l': l, 
            'kappa_xx': kappa_xx, 
            'kappa_yy': kappa_yy,
            'kappa_zz': kappa_zz,
            'kappa_perp': kappa_perp,
            'kappa_para': kappa_zz
        }
        return pd.DataFrame(d)
