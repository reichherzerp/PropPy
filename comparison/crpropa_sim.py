import crpropa as crp
import numpy as np
import pandas as pd


class CRPropa:
    def __init__(self, step_size = 10**11, traj_max = 10**14):
        # all simulation parameters
        self.energy = 3*10**16*crp.eV
        self.n_obs = 100
        self.n_particles = 10**3
        self.brms = crp.gauss
        self.l_max = 5*10**11 # [m]
        self.l_min = 5*10**9 # [m]
        self.n_wavemodes = 250
        self.step_size = step_size
        self.traj_max = traj_max
        self.file_name = 'data/raw_data/sim_result_'


    def set_energy(self, energy):
        self.energy = energy

    def set_n_obs(self, n_obs):
        self.n_obs = n_obs

    def set_n_particles(self, n_particles):
        self.n_particles = n_particles

    def set_brms(self, brms):
        self.brms = brms

    def set_l_max(self, l_max):
        self.l_max = l_max

    def set_l_min(self, l_min):
        self.l_min = l_min

    def set_n_wavemodes(self, n_wavemodes):
        self.n_wavemodes = n_wavemodes

    def set_step_size(self, step_size):
        self.step_size = step_size

    def set_traj_max(self, traj_max):
        self.traj_max = traj_max

    def set_file_name(self, file_name):
        self.file_name = file_name

    def sim(self):
        sim = crp.ModuleList()

        # point source settings
        source = crp.Source()
        source.add(crp.SourcePosition(crp.Vector3d(0)))
        source.add(crp.SourceParticleType(crp.nucleusId(1, 1)))
        source.add(crp.SourceEnergy(self.energy))
        source.add(crp.SourceIsotropicEmission())

        # magnetic field 
        b_field = crp.MagneticFieldList()
        turbulence_spectrum = crp.SimpleTurbulenceSpectrum(self.brms, self.l_min, self.l_max)
        turbulence = crp.PlaneWaveTurbulence(turbulence_spectrum, self.n_wavemodes)
        b_field.addField(turbulence)
        
        # propagation
        prop_bp = crp.PropagationBP(b_field, self.step_size)
        sim.add(prop_bp)
        maxTra = crp.MaximumTrajectoryLength(self.traj_max)
        sim.add(maxTra)

        # output
        output = crp.TextOutput(self.file_name+str(self.step_size/10**11)+'.txt', crp.Output.Trajectory3D)
        output.enable(output.SerialNumberColumn)

        # observer
        obs = crp.Observer()
        log = True
        obs.add(crp.ObserverTimeEvolution(self.step_size, self.traj_max, self.n_obs, log))
        obs.setDeactivateOnDetection(False)
        obs.onDetection(output)
        sim.add(obs)

        # run simulation     
        r_g = (self.energy/(crp.eV*crp.c_light*self.brms))
        l_c = self.l_max/5
        print('step_size/r_g = ', self.step_size/r_g)
        print('step_size/l_c = ', self.step_size/l_c)
        if self.step_size > r_g:
            print("warning: step size doesn't resolve gyromotion")
        if self.step_size > l_c:
            print("warning: step size doesn't resolve correlation length of turbulence")
        sim.setShowProgress(True)
        sim.run(source, self.n_particles, True)


    def load_data(self, x):
        dataI = pd.read_csv(x, names=['D', 'SN', 'ID', 'E', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz', 'SN0', 'SN1'], delimiter='\t', comment='#', usecols=["D", "X", "Y", "Z", "SN"])

        ### Convert data from Mpc to meters                                            
        dataI.X = dataI.X * 10**6*crp.pc
        dataI.Y = dataI.Y * 10**6*crp.pc
        dataI.Z = dataI.Z * 10**6*crp.pc
        dataI.D = dataI.D * 10**6*crp.pc
        ### Calculate Diffusion Coefficients
        dataI['X2D'] = dataI.X**2 / (dataI.D) * crp.c_light / 2.
        dataI['Y2D'] = dataI.Y**2 / (dataI.D) * crp.c_light / 2.
        dataI['Z2D'] = dataI.Z**2 / (dataI.D) * crp.c_light / 2.
        dataI['R2D'] = (dataI.X**2+dataI.Y**2+dataI.Z**2)**0.5
        ### Distance                                                            
        dataI.D = dataI.D
        print(len(dataI.D.values.tolist()))
        
        dataI = dataI.sort_values('D')
        return dataI


    def diffusion_coefficient_isotropic(self, data):
        kappa = []
        trajectory_lengths = self.return_l(data)
        for l in trajectory_lengths:
            dataI = data[data['D'] == l]
            kappa.append(np.mean(dataI.X2D.values + dataI.Y2D.values + dataI.Y2D.values)/3.0)
        return kappa


    def return_l(self, dataI):
        l = list(set(dataI.D.values.tolist()))
        return sorted(l)


    def analyze(self, step_size, file_name_output):
        data = self.load_data(self.file_name+str(step_size/10**11)+'.txt')
        kappa = self.diffusion_coefficient_isotropic(data)
        l = self.return_l(data)
        np.save(file_name_output+str(step_size/10**11)+'_l', np.array(l))
        np.save(file_name_output+str(step_size/10**11)+'_kappa.npy', np.array(kappa))