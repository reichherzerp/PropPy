"""Scripts to visualize simulated particle trajectories.

The simulated particles can be visualized here by showing their trajectories
or other properties.


    Typical usage example:

    import proppy as pp

    df = pd.read_pickle("data.pkl")
    df_time_evolution_observer = df.loc[df['radius'] == -1.0]
    tra = pp.Trajectory(df_time_evolution_observer)
    particle_ids = tra.get_particle_ids()

    nr_steps = 2*10**2
    tra.plot_trajectory('x', 'y', 'd', particle_ids[0], nr_steps, None)
    tra.plot_trajectory('d', 'z', 'd', particle_ids[0], nr_steps, None)
"""


import numpy as np
import matplotlib.pyplot as plt


class Trajectory():
    """Trajectory class for visualizing particle trajectories.

    After loading the simulation data, the simulated particles can be visualized 
    here by showing their trajectories or other properties.

    Attributes:
        dimensions: An int for defining the dimensions.
        df: A pandas dataftrame with the simulation data.
    """

    def __init__(self, df):
        print('init trajectory plotting class')
        self.df = df
        self.dimensions = 3


    def set_dimensions(self, dimensions):
        self.dimensions = dimensions

        
    def get_particle_ids(self):
        return list(map(int, (set(self.df['id']))))
        
        
    def plot_trajectory(self, x, y, c, particle_ids, nr_steps, file_name):
        """Plot the particle trajectory.
        
        Args:
            x: A String that specifies which column to plot on the x-axis.
            y: A String that specifies which column to plot on the y-axis.
            c: A String that specifies which column to use for color-coding the plot and the color-bar.
            particle_ids: List of particle ids that should be plotted.
            nr_steps: An int tha defines the number of steps to be plotted.
            file_name: String or None. If not None, the plot will be saved with the given String as a name.
        """
        if isinstance(particle_ids, int):
            # create a list if only one id was passed to particle_ids
            particle_ids = [particle_ids]
        plt.figure(figsize=(4.5,4))
        for particle_id in particle_ids:
            # filter the pandas dataframe for data of the current particle_id
            df_ids = self.df[self.df['id'] == particle_id]
            # filter the pandas dataframe for data of after each substep 
            # -> remove points from intermediate substeps
            df_ids = df_ids[df_ids['sub_step'] == self.dimensions-1]
            # color code the trajectory according to the distance travelled by the particle
            plt.scatter(df_ids[x][:nr_steps], df_ids[y][:nr_steps], s = 4, c = df_ids[c][:nr_steps], cmap = 'viridis')
        # plot colorbar next to plot
        cbar = plt.colorbar()
        cbar.set_label(c + ' [m]')
        plt.tight_layout()
        xlabel = self.get_label(x)
        ylabel = self.get_label(y)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        if file_name is not None:
            plt.savefig(file_name)
        plt.show()


    def get_label(self, axis):
        """If the axis shows a distance, add the unit [m] to the label.
        """
        label = axis
        if axis == 'x' or axis == 'y' or axis == 'z' or axis == 'd':
            label = label + ' [m]'
        return label
        

    def plot_trjectory_substeps(self, substep_0, substep_1, particle_ids, number_steps, file_name):
        """Function to visualize substeps. Two substeps can be passed, 
        which are visualized in the x-y plane. The lines from substep_0 
        to substep_1 as well as the lines substep_1 to aubstep_0 are marked in color
        """
        if isinstance(particle_ids, int):
            # create a list if only one id was passed to particle_ids
            particle_ids = [particle_ids]
        n = number_steps*self.dimensions+1 # number of substeps
        for particle_id in particle_ids:
            # filter the pandas dataframe for data of the current particle_id
            df_ids = self.df[self.df['id'] == particle_id][:n]
            df_substep_0 = df_ids[df_ids['sub_step'] == substep_0]
            df_substep_1 = df_ids[df_ids['sub_step'] == substep_1]
            df_steps = df_ids[df_ids['sub_step'] == self.dimensions-1]
            if len(df_substep_0) == 0 or len(df_substep_1) == 0:
                print('Error: Data has no substeps. Please observe substeps in simulation. Afterwards you can visualize them.')
                return
            
            x_substep_0 = df_substep_0['x']
            y_substep_0 = df_substep_0['y']
            x_substep_1 = df_substep_1['x']
            y_substep_1 = df_substep_1['y']
            # plot complete trajectory by connecting points after all substeps
            plt.figure(figsize=(4,4))
            plt.plot(df_steps['x'], df_steps['y'], c='grey', label='complete trajectory')
            # plot positions after substeps
            plt.scatter(x_substep_0[:-1], y_substep_0[:-1], c='dodgerblue', marker='^', s=15, label='after substep '+str(substep_0))
            plt.scatter(x_substep_1, y_substep_1, c='peachpuff', marker='d', s=15, label='after substep '+str(substep_1))
            # plot substep lines
            # generate the individual partial sections of the substeps between the points.
            # The origin gets the step number s = dimensions - 2. In 3 d the 0th step would be s = 2. 
            # Thus there are possibly more partial steps here than for 0 and 1, 
            # which is why this must be handled as follows:
            if len(x_substep_0) + 1 == len(x_substep_1):
                x_0_1 = np.vstack([x_substep_0[:],x_substep_1[:-1]])
                y_0_1 = np.vstack([y_substep_0[:],y_substep_1[:-1]])
                x_1_0 = np.vstack([x_substep_1[1:],x_substep_0[:]])
                y_1_0 = np.vstack([y_substep_1[1:],y_substep_0[:]])
            else:
                x_0_1 = np.vstack([x_substep_0[1:-1],x_substep_1[:-1]])
                y_0_1 = np.vstack([y_substep_0[1:-1],y_substep_1[:-1]])
                x_1_0 = np.vstack([x_substep_1[:],x_substep_0[:-1]])
                y_1_0 = np.vstack([y_substep_1[:],y_substep_0[:-1]])
            # legend for substeps
            plt.plot([x_substep_0.tolist()[0]], [y_substep_0.tolist()[0]], 'dodgerblue', ls=':', label='move '+str(substep_0)+'???'+str(substep_1))
            plt.plot([x_substep_1.tolist()[0]], [y_substep_1.tolist()[0]], 'peachpuff', ls='--', label='move '+str(substep_1)+'???'+str(substep_0))
            plt.plot(x_0_1, y_0_1, 'dodgerblue', ls=':')
            plt.plot(x_1_0, y_1_0, 'peachpuff', ls='--')
            
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.legend()
        plt.axis('square')
        plt.tight_layout()
        if file_name is not None:
            plt.savefig(file_name)
        plt.show()