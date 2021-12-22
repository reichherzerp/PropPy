import pandas as pd
import numpy as np
import math
pc = 3.086*10**16
Mpc = 10**6*pc
c = 2.99792*10**8
N_time = 100
number_particles = 10**3
N = number_particles

def load_data(x):
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
    ### Number of Gyrations                                                            
    dataI.D = dataI.D #/ (2 * math.pi * r)
    print(len(dataI.D.values.tolist()))
    
    dataI = dataI.sort_values('D')
    return dataI


def diffusion_coefficients(data):
    kappa_xx = []
    kappa_yy = []
    kappa_zz = []
    kappa_rr = []
    L = return_L(data)
    for l in L:
        dataI = data[data['D'] == l]
        kappa_xx.append(np.mean(dataI.X2D.values + dataI.Y2D.values)/2.0)
        kappa_zz.append(np.mean(dataI.Z2D.values))
    return kappa_xx, kappa_zz

def return_L(dataI):
    L = list(set(dataI.D.values.tolist()))
    return sorted(L)


def analyze_agn(file_name, file_name_output, diffusive):
    kappa_xx_averaged = np.zeros(N_time)
    kappa_zz_averaged = np.zeros(N_time)
    seeds = range(1)


    for seed in seeds:
        if diffusive:
            file_name = file_name + '_SED'
        dataLin = load_data(file_name+'.txt')
        kappa_xx, kappa_zz = diffusion_coefficients(dataLin)
        for i in range(N_time):
            kappa_xx_averaged[i] = kappa_xx_averaged[i] + kappa_xx[i]/len(seeds)
            kappa_zz_averaged[i] = kappa_zz_averaged[i] + kappa_zz[i]/len(seeds)

    L = return_L(dataLin)
    np.save(file_name_output+'_d', np.array(L))
    np.save(file_name_output+'_kappa_perp.npy', np.array(kappa_xx_averaged))
    np.save(file_name_output+'_kappa_para.npy', np.array(kappa_zz_averaged))

analyze_agn('sim_result', 'sim_result_ana', False)