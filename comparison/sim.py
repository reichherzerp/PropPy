import crpropa as crp
import numpy as np


sim = crp.ModuleList()

# source settings
# A point source at the origin is isotropically injecting protons.
source = crp.Source()
source.add(crp.SourcePosition(crp.Vector3d(0)))
source.add(crp.SourceParticleType(crp.nucleusId(1, 1)))
source.add(crp.SourceEnergy(3*10**16*crp.eV))
source.add(crp.SourceIsotropicEmission())

# magnetic field                                                                                              
b_field = crp.MagneticFieldList()
brms = 3*10**6*crp.muG
b_background = 10**7*crp.muG
l_max = 5*10**11 # [m]
l_min = 5*10**9 # [m]
turbulence_spectrum = crp.SimpleTurbulenceSpectrum(brms, l_min, l_max)
turbulence = crp.PlaneWaveTurbulence(turbulence_spectrum, Nm = 400)
b_field.addField(turbulence)
regular_field = crp.UniformMagneticField(crp.Vector3d(0,0,1)*b_background)
b_field.addField(regular_field)

# propagation
step_size = 10**10 #m
prop_bp = crp.PropagationBP(b_field, step_size)
sim.add(prop_bp)
t_max = 10**5*step_size
maxTra = crp.MaximumTrajectoryLength(t_max)
sim.add(maxTra)
