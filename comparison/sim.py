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
