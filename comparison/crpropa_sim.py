import crpropa as crp
import numpy as np
import pandas as pd
pc = 3.086*10**16
Mpc = 10**6*pc
c = 2.99792*10**8
N_time = 100
number_particles = 10**3
N = number_particles
file_name = 'data/sim_result_'

class CRPropa:
    def __init__(self, step_size = 10**11, traj_max = 10**16):
        # all simulation parameters
        self.energy = 3*10**16*crp.eV
        self.n_obs = 100
        self.brms = 10**6*crp.muG
        self.l_max = 5*10**11 # [m]
        self.l_min = 5*10**9 # [m]
        self.step_size = step_size
        self.traj_max = traj_max


    def crpropa_sim(self):
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
        turbulence = crp.PlaneWaveTurbulence(turbulence_spectrum, Nm = 100)
        b_field.addField(turbulence)
        
        # propagation
        prop_bp = crp.PropagationBP(b_field, self.step_size)
        sim.add(prop_bp)
        maxTra = crp.MaximumTrajectoryLength(self.traj_max)
        sim.add(maxTra)

        # output
        output = crp.TextOutput(file_name+str(self.step_size/10**11)+'.txt', crp.Output.Trajectory3D)
        output.enable(output.SerialNumberColumn)

        # observer
        obs = crp.Observer()
        obs.add(crp.ObserverTimeEvolution(self.step_size, self.traj_max, self.n_obs, log=1))
        obs.setDeactivateOnDetection(False)
        obs.onDetection(output)
        sim.add(obs)

        # run simulation                                               
        sim.setShowProgress(True)
        sim.run(source, number_particles, True)


    def load_data(self, x):
        dataI = pd.read_csv(x, names=['D', 'SN', 'ID', 'E', 'X', 'Y', 'Z', 'Px', 'Py', 'Pz', 'SN0', 'P0x', 'P0y', 'P0z', 'SN1'], delimiter='\t', comment='#', usecols=["D", "X", "Y", "Z", "\
    SN"])

        ### Convert data from Mpc to meters                                            
        dataI.X = dataI.X * Mpc
        dataI.Y = dataI.Y * Mpc
        dataI.Z = dataI.Z * Mpc
        dataI.D = dataI.D * Mpc
        ### Calculate Diffusion Coefficients
        dataI['X2D'] = dataI.X**2 / (dataI.D) * c / 2.
        dataI['Y2D'] = dataI.Y**2 / (dataI.D) * c / 2.
        dataI['Z2D'] = dataI.Z**2 / (dataI.D) * c / 2.
        dataI['R2D'] = (dataI.X**2+dataI.Y**2+dataI.Z**2)**0.5
        ### Distance                                                            
        dataI.D = dataI.D
        print(len(dataI.D.values.tolist()))
        
        dataI = dataI.sort_values('D')
        return dataI


    def diffusion_coefficients(self, data):
        kappa_xx = []
        kappa_yy = []
        kappa_zz = []
        kappa_rr = []
        L = self.return_L(data)
        for l in L:
            dataI = data[data['D'] == l]
            kappa_xx.append(np.mean(dataI.X2D.values + dataI.Y2D.values)/2.0)
            kappa_zz.append(np.mean(dataI.Z2D.values))
        return kappa_xx, kappa_zz

    def return_L(self, dataI):
        L = list(set(dataI.D.values.tolist()))
        return sorted(L)


    def analyze_agn(self, step_size, file_name_output):
        dataLin = self.load_data(file_name+str(step_size/10**11)+'.txt')
        kappa_perp, kappa_para = self.diffusion_coefficients(dataLin)
        L = self.return_L(dataLin)
        np.save(file_name_output+str(step_size/10**11)+'_d', np.array(L))
        np.save(file_name_output+str(step_size/10**11)+'_kappa_perp.npy', np.array(kappa_perp))
        np.save(file_name_output+str(step_size/10**11)+'_kappa_para.npy', np.array(kappa_para))