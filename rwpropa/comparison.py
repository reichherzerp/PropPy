import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
import pandas as pd
from scipy.optimize import curve_fit

class Comparison():

    def __init__(self, kappa_theory, lambda_theory, step_sizes, l_c, r_g, path_data_raw, path_data, path_figs, proppy_unit = 'm'):
        self.kappa_theory = kappa_theory
        self.lambda_theory = lambda_theory
        self.proppy_unit = proppy_unit
        self.step_sizes = step_sizes
        self.l_c = l_c
        self.r_g = r_g
        self.path_data = path_data
        self.path_data_raw = path_data_raw
        self.path_figs = path_figs
        self.load_sim_data()

    def load_sim_data(self):
        try:
            self.df_proppy_results = pd.read_pickle(self.path_data+'/proppy_sim_data.pkl')
        except:
            print("couldn't loade PropPy data")
        try:
            self.df_crp_ck_results = pd.read_pickle(self.path_data+'/crp_sim_data_CK_PW.pkl')
        except:
            print("couldn't loade CK data")
        try:
            self.df_crp_bp_pw_results = pd.read_pickle(self.path_data+'/crp_sim_data_BP_PW.pkl')
        except:
            print("couldn't loade BP PW data")
        try:
            self.df_crp_bp_grid_results = pd.read_pickle(self.path_data+'/crp_sim_data_BP_grid.pkl')
        except:
            print("couldn't loade BP grid data")
        try:
            self.df_crp_sde_results = pd.read_pickle(self.path_data+'/crp_sim_data_SDE_.pkl')
        except:
            print("couldn't loade SDE data")

        try:
            self.proppy_times = np.array(self.df_proppy_results['time'].values.tolist())
            if self.proppy_unit == 'm':
                self.proppy_step_sizes = self.df_proppy_results['step_size'] # [m]
                self.proppy_kappas = self.df_proppy_results['kappa'] # [m^2/s]
            if self.proppy_unit == 'km':
                self.proppy_step_sizes = self.df_proppy_results['step_size']*10**3 # [m]
                self.proppy_kappas = self.df_proppy_results['kappa']*10**6 # [m^2/s]
            if self.proppy_unit == 'pc':
                self.proppy_step_sizes = self.df_proppy_results['step_size']*3.086*10**16 # [m]
                self.proppy_kappas = self.df_proppy_results['kappa']*(3.086*10**16)**2 # [m^2/s]
        except:
            print('no PropPy data')
            self.proppy_times = np.array([])
            self.proppy_step_sizes = np.array([])
            self.proppy_kappas = np.array([])
            
        try:
            self.ck_times = np.array(self.df_crp_ck_results['time'].values.tolist())
            self.ck_step_sizes = self.df_crp_ck_results['step_size']
            self.ck_kappas = self.df_crp_ck_results['kappa']
        except:
            print('no ck data')
            self.ck_times = np.array([])
            self.ck_step_sizes = np.array([])
            self.ck_kappas = np.array([])
        try:
            self.bp_pw_times = np.array(self.df_crp_bp_pw_results['time'].values.tolist())
            self.bp_pw_step_sizes = self.df_crp_bp_pw_results['step_size']
            self.bp_pw_kappas = self.df_crp_bp_pw_results['kappa']
        except:
            print('no bp pw data')
            self.bp_pw_times = np.array([])
            self.bp_pw_step_sizes = np.array([])
            self.bp_pw_kappas = np.array([])
        try:
            self.bp_grid_times = np.array(self.df_crp_bp_grid_results['time'].values.tolist())
            self.bp_grid_step_sizes = self.df_crp_bp_grid_results['step_size']
            self.bp_grid_kappas = self.df_crp_bp_grid_results['kappa']
        except:
            print('no bp grid data')
            self.bp_grid_times = np.array([])
            self.bp_grid_step_sizes = np.array([])
            self.bp_grid_kappas = np.array([])
        try:
            self.sde_times = np.array(self.df_crp_sde_results['time'].values.tolist())
            self.sde_step_sizes = self.df_crp_sde_results['step_size']
            self.sde_kappas = self.df_crp_sde_results['kappa']
        except:
            print('no SDE data')
            self.sde_times = np.array([])
            self.sde_step_sizes = np.array([])
            self.sde_kappas = np.array([])


    def func_kappa_run(self, x, kappa, tau):
        return (1 - np.exp(-x/tau))*kappa


    def plot_running_diffusion_coefficients(self, d_theory=[1e17, 4e17], pp_step_max=0, times=[], kappas_ref=[], show_lambda=True):
        fig = plt.figure(figsize=(5,5))
        # set height ratios for subplots
        gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1]) 
        # remove vertical gap between subplots
        plt.subplots_adjust(hspace=.035)

        ### 1. subplot
        ax0 = plt.subplot(gs[0])

        plt.plot(d_theory, [self.kappa_theory*10**4,self.kappa_theory*10**4], color='k', linestyle=(0, (3, 1, 1, 1)), label='theory', zorder=-1)
        if show_lambda:
            plt.axvline(x=self.lambda_theory, color='k', linestyle='--', label='$\lambda_\mathrm{theory}$', zorder=-1)
        steps_proppy = []
        steps_ck = []
        steps_bp_pw = []
        steps_bp_grid = []
        steps_sde = []
        kappas_proppy = []
        kappas_ck = []
        kappas_bp_pw = []
        kappas_bp_grid = []
        kappas_sde = []

        plt.xlim([0.5*self.step_sizes[-1], d_theory[-1]*1.5])
        for i, step_size in enumerate(self.step_sizes):
            color = plt.cm.viridis(np.linspace(0, 1, len(self.step_sizes))[i])
            n_max = -1
            try:
                proppy_l = np.load(self.path_data+'/sim_result_proppy_stepsize_'+str(step_size/10**11)+'_l.npy')
                proppy_kappa = np.load(self.path_data+'/sim_result_proppy_stepsize_'+str(step_size/10**11)+'_kappa.npy')
                if self.proppy_unit == 'm':
                    pass
                elif self.proppy_unit == 'km':
                    proppy_l = proppy_l*10**3 # [m]
                    proppy_kappa = proppy_kappa*10**6 # [m^2/s]
                elif self.proppy_unit == 'pc':
                    proppy_l = proppy_l*3.086*10**16 # [m]
                    proppy_kappa = proppy_kappa*(3.086*10**16)**2 # [m^2/s]
                if pp_step_max == 0 or step_size < self.lambda_theory:
                    ax0.plot(proppy_l[:n_max], np.array(proppy_kappa[:n_max])*10**4, color='red', ls='-', zorder=2, lw=2) 
                steps_proppy.append(step_size)
                kappas_proppy.append(np.mean(proppy_kappa[-10:]))
            except Exception as e: 
                print(e)
                print('no data for PropPy')
            
            try:
                crp_l = np.load(self.path_data+'/sim_result_crp_BP_PW_stepsize_'+str(step_size/10**11)+'_l.npy')
                crp_kappa = np.load(self.path_data+'/sim_result_crp_BP_PW_stepsize_'+str(step_size/10**11)+'_kappa.npy')
                ax0.plot(crp_l[:n_max], np.array(crp_kappa[:n_max])*10**4, color=color, ls=(0, (1, 1)), lw=2, zorder=4)
                steps_bp_pw.append(step_size)
                kappas_bp_pw.append(np.mean(crp_kappa[-10:]))
            except Exception as e: 
                print(e)
                print('no data for BP pw')
            try:
                crp_l = np.load(self.path_data+'/sim_result_crp_BP_grid_stepsize_'+str(step_size/10**11)+'_l.npy')
                crp_kappa = np.load(self.path_data+'/sim_result_crp_BP_grid_stepsize_'+str(step_size/10**11)+'_kappa.npy')
                ax0.plot(crp_l[:n_max], np.array(crp_kappa[:n_max])*10**4, color=color, ls=(0, (1, 8)), lw=2, zorder=4)
                steps_bp_grid.append(step_size)
                kappas_bp_grid.append(np.mean(crp_kappa[-10:]))
            except Exception as e: 
                print(e)
                print('no data for BP grid')
            
            try:
                crp_l = np.load(self.path_data+'/sim_result_crp_CK_PW_stepsize_'+str(step_size/10**11)+'_l.npy')
                crp_kappa = np.load(self.path_data+'/sim_result_crp_CK_PW_stepsize_'+str(step_size/10**11)+'_kappa.npy')
                ax0.plot(crp_l[:n_max], np.array(crp_kappa[:n_max])*10**4, color=color, ls='-.', lw=2, zorder=4)
                steps_ck.append(step_size)
                kappas_ck.append(np.mean(crp_kappa[-10:]))
            except Exception as e: 
                print(e)
                print('no data for CK')

            try:
                crp_l = np.load(self.path_data+'/sim_result_crp_SDE__stepsize_'+str(step_size/10**11)+'_l.npy')
                crp_kappa = np.load(self.path_data+'/sim_result_crp_SDE__stepsize_'+str(step_size/10**11)+'_kappa.npy')
                ax0.plot(crp_l[:n_max], np.array(crp_kappa[:n_max])*10**4, color=color, ls='--', lw=2, zorder=4)
                steps_sde.append(step_size)
                kappas_sde.append(np.mean(crp_kappa[-10:]))
            except Exception as e: 
                print(e)
                print('no data for SDE')

        # colorbar
        plt.scatter(np.zeros(len(self.step_sizes)), np.zeros(len(self.step_sizes)), c=self.step_sizes, cmap='viridis', norm=matplotlib.colors.LogNorm())
        plt.colorbar(label='step sizes [m]', pad=0.01)

        #legend
        plt.plot([0,0], [0,0], c='grey', ls=':', label='CRPropa (BP) [PW]', lw=2)
        plt.plot([0,0], [0,0], c='grey', ls=(0, (1, 8)), label='CRPropa (BP) [grid]', lw=2)
        plt.plot([0,0], [0,0], c='grey', ls='-.', label='CRPropa (CK) [PW]', lw=2)
        plt.plot([0,0], [0,0], c='grey', ls='--', label='CRPropa (SDE)', lw=2)
        plt.plot([0,0], [0,0], c='red', ls='-', label='PropPy', lw=2)

        ax0.tick_params(labelbottom=False)
        ax0.tick_params(top=True, right=True, direction='in', which='both')
        ax0.loglog()
        ax0.set_ylabel('running $\kappa$ [cm$^2$/s]')
        plt.legend()


        ax1 = plt.subplot(gs[1], sharex = ax0)
        ax1.set_xlabel('trajectory length [m]')
        ax1.set_ylabel('deviation [%]')

        if len(times) > 0:
            # Running diffusion coefficients can be modelled via
            # \begin{equation}
            # \kappa (t) = \kappa \cdot \left(1 - \mathrm{e}^{-t/\tau}\right),
            # \end{equation}
            # where $\kappa$ is the converged diffusion coefficient and $\tau$ can be deduced from the initial ballistic phase via
            # \begin{equation}
            # \tau = \frac{t}{\kappa(t)}
            # \end{equation}
            # for small $t$.
            # Note that this expression was derived by using the Taylor series of $\kappa (t)$ for small times
            # \begin{equation}
            # \kappa (t) \approx \frac{t}{\tau}.
            # \end{equation}
            tau = times[0] / kappas_ref[0] * kappas_ref[-1]
        

        for i, step_size in enumerate(self.step_sizes[:-2]):
            color = plt.cm.viridis(np.linspace(0, 1, len(self.step_sizes))[i])
            n_max = -1
            try:
                proppy_l = np.load(self.path_data+'/sim_result_proppy_stepsize_'+str(step_size/10**11)+'_l.npy')
                proppy_kappa = np.load(self.path_data+'/sim_result_proppy_stepsize_'+str(step_size/10**11)+'_kappa.npy')
                if self.proppy_unit == 'm':
                    pass
                elif self.proppy_unit == 'km':
                    proppy_l = proppy_l*10**3 # [m]
                    proppy_kappa = proppy_kappa*10**6 # [m^2/s]
                elif self.proppy_unit == 'pc':
                    proppy_l = proppy_l*3.086*10**16 # [m]
                    proppy_kappa = proppy_kappa*(3.086*10**16)**2 # [m^2/s]

                if pp_step_max == 0 or step_size < self.lambda_theory:
                    kappa_theory = self.func_kappa_run(proppy_l, kappas_ref[-1], tau)
                    kappa_dev = np.abs(kappa_theory - proppy_kappa) / kappa_theory * 100
                    plt.plot(proppy_l, kappa_dev, color='red', ls='-', zorder=2, lw=2)
            except:
                print('no data for PropPy')


            try:
                crp_l = np.load(self.path_data+'/sim_result_crp_BP_PW_stepsize_'+str(step_size/10**11)+'_l.npy')
                crp_kappa = np.load(self.path_data+'/sim_result_crp_BP_PW_stepsize_'+str(step_size/10**11)+'_kappa.npy')
                kappa_theory = self.func_kappa_run(crp_l, kappas_ref[-1], tau)
                kappa_dev = np.abs(kappa_theory - crp_kappa) / kappa_theory * 100
                plt.plot(crp_l, kappa_dev, color=color, ls=(0, (1, 1)), lw=2)
            except:
                print('no data for BP pw')
            

            try:
                crp_l = np.load(self.path_data+'/sim_result_crp_SDE__stepsize_'+str(step_size/10**11)+'_l.npy')
                crp_kappa = np.load(self.path_data+'/sim_result_crp_SDE__stepsize_'+str(step_size/10**11)+'_kappa.npy')
                kappa_theory = self.func_kappa_run(crp_l, kappas_ref[-1], tau)
                kappa_dev = np.abs(kappa_theory - crp_kappa) / kappa_theory * 100
                plt.plot(crp_l, kappa_dev, color=color, ls='--', lw=2)
            except:
                print('no data for SDE')


        plt.yscale('log')

        self.color_bar_invisible()
        ax1.tick_params(top=True, right=True, direction='in', which='both')
        plt.axhline(y=1, color='grey', ls='-', zorder=-4)
        
        plt.savefig(self.path_figs+'/running_kappa.pdf', bbox_inches='tight', pad_inches=0.02)
        plt.show()


        # second figure

        fig, ax1 = plt.subplots(figsize=(5,3.5))
        plt.scatter(steps_proppy, kappas_proppy, label='PropPy', marker='s', color='green')
        plt.scatter(steps_ck, kappas_ck, label='CRPropa (CK) [PW]', color='r')
        plt.scatter(steps_bp_pw, kappas_bp_pw, label='CRPropa (BP) [PW]', marker='d', color='k')
        plt.scatter(steps_bp_grid, kappas_bp_grid, label='CRPropa (BP) [grid]', marker='*', color='purple')
        plt.scatter(steps_sde, kappas_sde, label='CRPropa (SDE)', marker='^', color='blue')
        plt.axvline(x=self.l_c, label='$l_\mathrm{c}$', color='grey', ls=':')
        plt.axvline(x=self.r_g*2*3.14, label='$2\pi\, r_\mathrm{g}$', color='grey', ls='--')
        plt.axvline(x=self.lambda_theory, label='$\lambda_\mathrm{theory}$', color='grey', ls='-.')
        plt.axhline(y=self.kappa_theory, color='k', linestyle='-', label='theory')
        plt.legend()
        plt.loglog()
        plt.show()


    def color_bar_invisible(self):
        # invisible colorbar for correct spacing
        cb = plt.colorbar(pad=0.01, aspect=1000)
        cb.set_label('', fontsize=10, labelpad=16)
        cb.ax.tick_params(size=0) #Remove ticks
        cb.outline.set_visible(False)
        cb.set_ticks([])


    def plot_kappa_convergence_tests(self, ylabel="$\kappa$ [m$^2$/s]", ylabel2="log deviation",  lambda_theory=True, kappa_mean_seeds=0, kappa_mean_seeds_err=0):
        fig = plt.figure(figsize=(5,5))
        # set height ratios for subplots
        gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1]) 
        # remove vertical gap between subplots
        plt.subplots_adjust(hspace=.035)

        ### 1. subplot
        ax0 = plt.subplot(gs[0])

        zs = np.concatenate([self.proppy_times, self.ck_times, self.bp_pw_times, self.bp_grid_times, self.sde_times], axis=0)
        min_, max_ = zs.min(), zs.max()
        plt.scatter(self.proppy_step_sizes, self.proppy_kappas, c=self.proppy_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='s')
        plt.clim(min_, max_)
        plt.scatter(self.ck_step_sizes, self.ck_kappas, c=self.ck_times, cmap='viridis', norm=matplotlib.colors.LogNorm())
        plt.clim(min_, max_)
        plt.scatter(self.bp_pw_step_sizes, self.bp_pw_kappas, c=self.bp_pw_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='d')
        plt.clim(min_, max_)
        plt.scatter(self.bp_grid_step_sizes, self.bp_grid_kappas, c=self.bp_grid_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='*')
        plt.clim(min_, max_)
        plt.scatter(self.sde_step_sizes, self.sde_kappas, c=self.sde_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='^')
        plt.clim(min_, max_)
        plt.colorbar(label='simulation time [s]', pad=0.01)
        plt.loglog()
        if lambda_theory:
            plt.axvline(x=self.lambda_theory, label='$\lambda_\mathrm{theory}$', color='grey', ls='-.', zorder=-1)
        plt.axvline(x=self.l_c, label='$l_\mathrm{c}$', color='grey', ls=':', zorder=-1)
        plt.axvline(x=self.r_g*2*3.14, label='$2\pi\, r_\mathrm{g}$', color='grey', ls='--', zorder=-1)
        plt.axhline(y=self.kappa_theory, color='k', linestyle='-', label='theory', zorder=-1)

        # legend
        plt.scatter([0],[0], label='PropPy', marker='s', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (CK) [PW]', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (BP) [PW]', marker='d', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (BP) [grid]', marker='*', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (SDE)', marker='^', color='grey', zorder=-1)

        plt.xlabel('step size [m]')
        plt.ylabel(ylabel)
        plt.legend()

        ax0.tick_params(labelbottom=False)
        ax0.tick_params(top=True, right=True, direction='in', which='both')

        ax1 = plt.subplot(gs[1], sharex = ax0)

        err_proppy = np.abs(np.log10(self.proppy_kappas)-np.log10(self.kappa_theory))
        err_crp_ck = np.abs(np.log10(self.ck_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_pw = np.abs(np.log10(self.bp_pw_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_grid = np.abs(np.log10(self.bp_grid_kappas)-np.log10(self.kappa_theory))
        err_crp_sde = np.abs(np.log10(self.sde_kappas)-np.log10(self.kappa_theory))
        
        plt.scatter(self.proppy_step_sizes, err_proppy, c=self.proppy_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='s')
        plt.clim(min_, max_)
        plt.scatter(self.ck_step_sizes, err_crp_ck, c=self.ck_times, cmap='viridis', norm=matplotlib.colors.LogNorm())
        plt.clim(min_, max_)
        plt.scatter(self.bp_pw_step_sizes, err_crp_bp_pw, c=self.bp_pw_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='d')
        plt.clim(min_, max_)
        plt.scatter(self.bp_grid_step_sizes, err_crp_bp_grid, c=self.bp_grid_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='*')
        plt.clim(min_, max_)
        plt.scatter(self.sde_step_sizes, err_crp_sde, c=self.sde_times, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='^')
        plt.clim(min_, max_)

        if lambda_theory:
            plt.axvline(x=self.lambda_theory, label='$\lambda_\mathrm{theory}$', color='grey', ls='-.', zorder=-1)
        plt.axvline(x=self.l_c, label='$l_\mathrm{c}$', color='grey', ls=':', zorder=-1)
        plt.axvline(x=self.r_g*2*3.14, label='$2\pi\, r_\mathrm{g}$', color='grey', ls='--', zorder=-1)

        plt.axhline(y=0, color='k', linestyle='-', zorder=-2)

        self.color_bar_invisible()
        ax1.tick_params(top=True, right=True, direction='in', which='both')

        plt.ylabel(ylabel2)

        plt.savefig(self.path_figs+'/kappa_vs_stepsize.pdf', bbox_inches='tight', pad_inches=0.02)
        plt.show()


    def plot_particle_distributions_crp(self, file_name, prop_type='BP', unit='m'):
        fig = plt.figure(figsize=(5,3.5))
        for i, step_size in enumerate(self.step_sizes):
            
            dataI = pd.read_csv(file_name+'_stepsize_'+str(step_size/10**11)+'.txt', names=['D', 'SN', 'ID', 'E', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz', 'SN0', 'SN1'], delimiter='\t', comment='#', usecols=["D", "X", "Y", "Z", "SN"])
            dataI = dataI.sort_values('D')
            max_l = max(dataI['D'].values.tolist())
            dataI = dataI[dataI['D'] == max_l]   
            color = plt.cm.viridis(np.linspace(0, 1, len(self.step_sizes))[i])
            xs = np.array(dataI.X.values.tolist())
            ys = np.array(dataI.Y.values.tolist())
            zs = np.array(dataI.Z.values.tolist())
            if unit == 'Mpc':
                scale = 3.086e22
            else:
                scale = 1
            plt.hist(np.concatenate((xs,ys,zs), axis=None)*scale, bins=50, range=[-1e17/5, 1e17/5], histtype=u'step', edgecolor=color, linewidth=1., facecolor="None")
            
                 
        # colorbar
        plt.scatter(np.zeros(len(self.step_sizes)), np.zeros(len(self.step_sizes)), s=0.0001, c=self.step_sizes, cmap='viridis', norm=matplotlib.colors.LogNorm())
        plt.colorbar(label='step sizes [m]')

        plt.xlabel('position $x_i$ [m]')
        plt.ylabel('# particles')
        plt.title('$l_\mathrm{traj} = 10^{17}$ m ('+prop_type+')')

        #plt.yscale('log')
        plt.savefig(self.path_figs+'/particle_distributions_'+prop_type+'.pdf', bbox_inches='tight', pad_inches=0.02) 
        plt.show()


    def plot_particle_distributions_proppy(self):
        prop_type = 'PropPy'
        unit = 'm'
        fig = plt.figure(figsize=(5,3.5))
        for i, step_size in enumerate(self.step_sizes):
            if i == len(self.step_sizes)-1:
                continue
            try:
                data = pd.read_pickle(self.path_data_raw+'/proppy_stepsize_'+str(step_size/10**11)+'.pkl')
                data = data.sort_values('d')
                max_l = max(data['d'].values.tolist())
                data = data[data['d'] == max_l]        
                bins = 10
                color = plt.cm.viridis(np.linspace(0, 1, len(self.step_sizes))[i])
                xs = np.array(data.x.values.tolist())
                ys = np.array(data.y.values.tolist())
                zs = np.array(data.z.values.tolist())
                if unit == 'Mpc':
                    scale = 3.086e22
                else:
                    scale = 1
                plt.hist(np.concatenate((xs,ys,zs), axis=None)*scale, bins=50, range=[-1e17/5, 1e17/5], histtype=u'step', edgecolor=color, linewidth=1., facecolor="None")
            except:
                pass
        # colorbar
        plt.scatter(np.zeros(len(self.step_sizes)), np.zeros(len(self.step_sizes)), c=self.step_sizes, cmap='viridis', s=0.00001, norm=matplotlib.colors.LogNorm())
        plt.colorbar(label='step sizes [m]')

        plt.xlabel('position $x_i$ [m]')
        plt.ylabel('# particles')
        plt.title('$l_\mathrm{traj} = 10^{17}$ m ('+prop_type+')')

        plt.savefig(self.path_figs+'/particle_distributions_'+prop_type+'.pdf', bbox_inches='tight', pad_inches=0.02) 
        plt.show()


    def plot_kappa_vs_time_steps(self):
        fig = plt.figure(figsize=(5,3.5))
        zs = np.concatenate([self.proppy_step_sizes, self.ck_step_sizes, self.bp_pw_step_sizes, self.bp_grid_step_sizes, self.sde_step_sizes], axis=0)
        min_, max_ = zs.min(), zs.max()
        plt.scatter(self.proppy_times, self.proppy_kappas, c=self.proppy_step_sizes, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='s')
        plt.clim(min_, max_)
        plt.scatter(self.ck_times, self.ck_kappas, c=self.ck_step_sizes, cmap='viridis', norm=matplotlib.colors.LogNorm())
        plt.clim(min_, max_)
        plt.scatter(self.bp_pw_times, self.bp_pw_kappas, c=self.bp_pw_step_sizes, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='d')
        plt.clim(min_, max_)
        plt.scatter(self.bp_grid_times, self.bp_grid_kappas, c=self.bp_grid_step_sizes, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='*')
        plt.clim(min_, max_)
        plt.scatter(self.sde_times, self.sde_kappas, c=self.sde_step_sizes, cmap='viridis', norm=matplotlib.colors.LogNorm(), marker='^')
        plt.clim(min_, max_)
        plt.colorbar(label='step size [s]')
        plt.loglog()
        plt.axhline(y=self.kappa_theory, color='k', linestyle='-', label='theory', zorder=-1)

        # legend
        plt.scatter([0],[0], label='PropPy', marker='s', color='grey')
        plt.scatter([0],[0], label='CRPropa (CK) [PW]', color='grey')
        plt.scatter([0],[0], label='CRPropa (BP) [PW]', marker='d', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (BP) [grid]', marker='*', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (SDE)', marker='^', color='grey')

        plt.xlabel('simulation time [s]')
        plt.ylabel('$\kappa$ [m$^2$/s]')
        plt.legend(loc = 'center left')
        plt.savefig(self.path_figs+'/kappa_vs_time_steps.pdf', bbox_inches='tight', pad_inches=0.02)
        plt.show()


    def plot_kappa_vs_time_deviation(self):
        fig = plt.figure(figsize=(5,3.5))
        err_proppy = np.abs(np.log10(self.proppy_kappas)-np.log10(self.kappa_theory))
        err_crp_ck = np.abs(np.log10(self.ck_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_pw = np.abs(np.log10(self.bp_pw_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_grid = np.abs(np.log10(self.bp_grid_kappas)-np.log10(self.kappa_theory))
        err_crp_sde = np.abs(np.log10(self.sde_kappas)-np.log10(self.kappa_theory))
        zs = np.concatenate([err_proppy, err_crp_ck, err_crp_bp_pw, err_crp_bp_grid, err_crp_sde], axis=0)
        zs = zs[~np.isnan(zs)]
        min_, max_ = zs.min(), zs.max()
        plt.scatter(self.proppy_times, self.proppy_kappas, c=err_proppy, cmap='viridis', marker='s')
        plt.clim(min_, max_)
        plt.scatter(self.ck_times, self.ck_kappas, c=err_crp_ck, cmap='viridis')
        plt.clim(min_, max_)
        plt.scatter(self.bp_pw_times, self.bp_pw_kappas, c=err_crp_bp_pw, cmap='viridis', marker='d')
        plt.clim(min_, max_)
        plt.scatter(self.bp_grid_times, self.bp_grid_kappas, c=err_crp_bp_grid, cmap='viridis', marker='*')
        plt.clim(min_, max_)
        plt.scatter(self.sde_times, self.sde_kappas, c=err_crp_sde, cmap='viridis', marker='^')
        plt.clim(min_, max_)
        plt.colorbar(label='deviation = |log($\kappa_\mathrm{sim}$) / log($\kappa_\mathrm{theory}$)|')
        plt.loglog()
        plt.axhline(y=self.kappa_theory, color='k', linestyle='-', label='theory', zorder=-1)

        # legend
        plt.scatter([0],[0], label='PropPy', marker='s', color='grey')
        plt.scatter([0],[0], label='CRPropa (CK) [PW]', color='grey')
        plt.scatter([0],[0], label='CRPropa (BP) [PW]', marker='d', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (BP) [grid]', marker='*', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (SDE)', marker='^', color='grey')

        plt.xlabel('simulation time [s]')
        plt.ylabel('$\kappa$ [m$^2$/s]')
        plt.legend(loc='center left')
        plt.savefig(self.path_figs+'/kappa_vs_time_deviation.pdf', bbox_inches='tight', pad_inches=0.02)
        plt.show()


    def plot_deviation_vs_time_steps(self, day=True):
        fig = plt.figure(figsize=(5,3.5))
        plt.axhline(y=0, color='k', linestyle='-', zorder=-1, label='theory')
        plt.axvline(x=1, color='grey', linestyle=(0, (5, 0.4)), zorder=-1, label='1 sec')
        plt.axvline(x=60, color='grey', linestyle=(0, (5, 2.5)), zorder=-1, label='1 min')
        plt.axvline(x=60*60, color='grey', linestyle=(0, (5, 7)), zorder=-1, label='1 hour')
        if day:
            plt.axvline(x=60*60*24, color='grey', linestyle=(0, (5, 13)), zorder=-1, label='1 day')
    
        err_proppy = np.abs(np.log10(self.proppy_kappas)-np.log10(self.kappa_theory))
        err_crp_ck = np.abs(np.log10(self.ck_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_pw = np.abs(np.log10(self.bp_pw_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_grid = np.abs(np.log10(self.bp_grid_kappas)-np.log10(self.kappa_theory))
        err_crp_sde = np.abs(np.log10(self.sde_kappas)-np.log10(self.kappa_theory)) 
        zs = np.concatenate([self.proppy_step_sizes, self.ck_step_sizes, self.bp_pw_step_sizes, self.bp_grid_step_sizes, self.sde_step_sizes], axis=0)
        min_, max_ = zs.min(), zs.max()
        plt.scatter(self.proppy_times, err_proppy, c=self.proppy_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='s')
        plt.clim(min_, max_)
        plt.scatter(self.ck_times, err_crp_ck, c=self.ck_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis')
        plt.clim(min_, max_)
        plt.scatter(self.bp_pw_times, err_crp_bp_pw, c=self.bp_pw_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='d')
        plt.clim(min_, max_)
        plt.scatter(self.bp_grid_times, err_crp_bp_grid, c=self.bp_grid_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='*')
        plt.clim(min_, max_)
        plt.scatter(self.sde_times, err_crp_sde, c=self.sde_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='^')
        plt.clim(min_, max_)
        plt.colorbar(label='step size [m]')
        plt.xscale('log')

        # legend
        plt.scatter([0],[0], label='PropPy', marker='s', color='grey')
        plt.scatter([0],[0], label='CRPropa (CK) [PW]', color='grey')
        plt.scatter([0],[0], label='CRPropa (BP) [PW]', marker='d', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (BP) [grid]', marker='*', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (SDE)', marker='^', color='grey')

        plt.xlabel('simulation time [s]')
        plt.ylabel('deviation = |log($\kappa_\mathrm{sim}$) / log($\kappa_\mathrm{theory}$)|')
        plt.legend(loc='upper right')
        plt.savefig(self.path_figs+'/deviation_vs_time_steps.pdf', bbox_inches='tight', pad_inches=0.02)
        plt.show()

    
    def plot_time_vs_deviation_steps(self, day=True):
        fig = plt.figure(figsize=(5,3.5))

        plt.axvline(x=0, color='k', linestyle='-', zorder=-1, label='theory')
        plt.axhline(y=1, color='grey', linestyle=(0, (5, 0.4)), zorder=-1, label='1 sec')
        plt.axhline(y=60, color='grey', linestyle=(0, (5, 2.5)), zorder=-1, label='1 min')
        plt.axhline(y=60*60, color='grey', linestyle=(0, (5, 7)), zorder=-1, label='1 hour')
        if day:
            plt.axhline(y=60*60*24, color='grey', linestyle=(0, (5, 13)), zorder=-1, label='1 day')

        err_proppy = np.abs(np.log10(self.proppy_kappas)-np.log10(self.kappa_theory))
        err_crp_ck = np.abs(np.log10(self.ck_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_pw = np.abs(np.log10(self.bp_pw_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_grid = np.abs(np.log10(self.bp_grid_kappas)-np.log10(self.kappa_theory))
        err_crp_sde = np.abs(np.log10(self.sde_kappas)-np.log10(self.kappa_theory)) 
        zs = np.concatenate([self.proppy_step_sizes, self.ck_step_sizes, self.bp_pw_step_sizes, self.bp_grid_step_sizes, self.sde_step_sizes], axis=0)
        min_, max_ = zs.min(), zs.max()
        plt.scatter(err_proppy, self.proppy_times, c=self.proppy_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='s')
        plt.clim(min_, max_)
        plt.scatter(err_crp_ck, self.ck_times, c=self.ck_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis')
        plt.clim(min_, max_)
        plt.scatter(err_crp_bp_pw, self.bp_pw_times, c=self.bp_pw_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='d')
        plt.clim(min_, max_)
        plt.scatter(err_crp_bp_grid, self.bp_grid_times, c=self.bp_grid_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='d')
        plt.clim(min_, max_)
        plt.scatter(err_crp_sde, self.sde_times, c=self.sde_step_sizes, norm=matplotlib.colors.LogNorm(), cmap='viridis', marker='^')
        plt.clim(min_, max_)
        plt.colorbar(label='step size [m]')
        plt.yscale('log')

        # legend
        plt.scatter([0],[0], label='PropPy', marker='s', color='grey')
        plt.scatter([0],[0], label='CRPropa (CK) [PW]', color='grey')
        plt.scatter([0],[0], label='CRPropa (BP) [PW]', marker='d', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (BP) [grid]', marker='*', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (SDE)', marker='^', color='grey')

        plt.ylabel('simulation time [s]')
        plt.xlabel('deviation = |log($\kappa_\mathrm{sim}$) / log($\kappa_\mathrm{theory}$)|')
        plt.legend(loc='upper right', ncol=2)
        plt.savefig(self.path_figs+'/time_vs_deviation_steps.pdf', bbox_inches='tight', pad_inches=0.02)
        plt.show()


    def plot_time_vs_steps_deviation(self, day = True, lambda_theory=True):
        fig = plt.figure(figsize=(5,3.5))
        err_proppy = np.abs(np.log10(self.proppy_kappas)-np.log10(self.kappa_theory))
        err_crp_ck = np.abs(np.log10(self.ck_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_pw = np.abs(np.log10(self.bp_pw_kappas)-np.log10(self.kappa_theory))
        err_crp_bp_grid = np.abs(np.log10(self.bp_grid_kappas)-np.log10(self.kappa_theory))
        err_crp_sde = np.abs(np.log10(self.sde_kappas)-np.log10(self.kappa_theory))
        zs = np.concatenate([err_proppy, err_crp_ck, err_crp_bp_pw, err_crp_bp_grid, err_crp_sde], axis=0)
        zs = zs[~np.isnan(zs)]
        min_, max_ = zs.min(), zs.max()
        plt.scatter(self.proppy_step_sizes, self.proppy_times, c=err_proppy, cmap='inferno', marker='s')
        plt.clim(min_, max_)
        plt.scatter(self.ck_step_sizes, self.ck_times, c=err_crp_ck, cmap='inferno')
        plt.clim(min_, max_)
        plt.scatter(self.bp_pw_step_sizes, self.bp_pw_times, c=err_crp_bp_pw, cmap='inferno', marker='d')
        plt.clim(min_, max_)
        plt.scatter(self.bp_grid_step_sizes, self.bp_grid_times, c=err_crp_bp_grid, cmap='inferno', marker='*')
        plt.clim(min_, max_)
        plt.scatter(self.sde_step_sizes, self.sde_times, c=err_crp_sde, cmap='inferno', marker='^')
        plt.clim(min_, max_)
        plt.colorbar(label='deviation = |log($\kappa_\mathrm{sim}$) / log($\kappa_\mathrm{theory}$)|')
        plt.loglog()
        
        # references for simulation time
        plt.axhline(y=1, color='grey', linestyle=(0, (5, 0.4)), zorder=-1, label='1 sec')
        plt.axhline(y=60, color='grey', linestyle=(0, (5, 2.5)), zorder=-1, label='1 min')
        plt.axhline(y=60*60, color='grey', linestyle=(0, (5, 7)), zorder=-1, label='1 hour')
        if day:
            plt.axhline(y=60*60*24, color='grey', linestyle=(0, (5, 13)), zorder=-1, label='1 day')

        # references for step sizes
        if lambda_theory:
            plt.axvline(x=self.lambda_theory, label='$\lambda_\mathrm{theory}$', color='grey', linestyle=(0, (5, 7)), zorder=-1)
        plt.axvline(x=self.l_c, label='$l_\mathrm{c}$', color='grey', linestyle=(0, (5, 0.4)), zorder=-1)
        plt.axvline(x=self.r_g*2*3.14, label='$2\pi\, r_\mathrm{g}$', color='grey', linestyle=(0, (5, 2.5)), zorder=-1)

        # legend
        plt.scatter([0],[0], label='PropPy', marker='s', color='grey')
        plt.scatter([0],[0], label='CRPropa (CK) [PW]', color='grey')
        plt.scatter([0],[0], label='CRPropa (BP) [PW]', marker='d', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (BP) [grid]', marker='*', color='grey', zorder=-1)
        plt.scatter([0],[0], label='CRPropa (SDE)', marker='^', color='grey')

        plt.xlabel('step sizes [m]')
        plt.ylabel('simulation time [s]')
        plt.legend(loc='center left')
        plt.savefig(self.path_figs+'/time_vs_steps_deviation.pdf', bbox_inches='tight', pad_inches=0.02)
        plt.show()