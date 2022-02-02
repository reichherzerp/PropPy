# PropPy

PropPy is a novel high-performance software tool to propagate individual charged particles in turbulent magnetic fields via a correlated random walk, while statistically satisfying the particle distribution described by the transport equation. This novel method is superior to diffusive propagation codes for compact objects, such as active galactic nuclei (AGN) blobs, pulsars, and supernovas because of its ability to describe the initial transport correctly. This novel method solves the problem of overestimating particle distributions in outer regions of the objects that exist in diffusion codes. Statistical properties of particles such as distributions, escape times, and running diffusion coefficients are in agreement with ballistic methods within the established propagation code CRPropa, while decreasing the computation time significantly.

# Getting started


### Installation
Simply clone the repository and import the package `import proppy as pp` while specifying the correct path.


### Example notebooks
- Start by having a look at the [basic simulation setup](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/tutorials/Tutorial_0---Template_Simulation.ipynb), where the important modules are introduced and basic simulation setup is shown. 

- Details of the individual modules are in [detailed simulation setup](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/tutorials/Tutorial_1---Detailed_Simulation_Setup.ipynb).

- If you are already familiar with diffusion coefficients, it is advisable to also try [anisotropic diffusion](https://gitlab.ruhr-uni-bochum.de/tp4/diffusion/rwpropa/-/blob/master/rwpropa/tutorials/Tutorial_2---Anisotropic_Diffusion.ipynb), which occurs due to an ordered magnetic field.

- To get a better feeling of how the code works, you should visualize [particle trajectories](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/tutorials/Tutorial_3---Particle_Trajectories.ipynb) with a dedicated notebook.

- [Spherical observers](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/tutorials/Tutorial_4---Tutorial_Spherical_Observer.ipynb) are a special type of observers that are useful in many astrophysical environments.

- A special usecase for sphereical observers is studied in more detail in the conteyt of [AGN plasmoids](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/tutorials/Tutorial_5---Tutorial_AGN_Plasmoid_Example.ipynb).

The software is [validated against CRPropa 3.1.7](https://github.com/CRPropa/CRPropa3/tree/3.1.7), which is an established code for charged particle propagation and interaction. In comparison, the strengths of RWPropa stand out, especially the good accuracy with the approach of solving the equation of motion. However, RWPropa is much faster.

Dedicated comparisons are shown in the folder [comparison](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/tree/master/comparison):
- [Active galactice nuclei jet plasmoids (ballistic case)](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/comparison/Comparison%20in%20compact%20astrophysical%20sources%201e14m%20trajectory.ipynb)
- [Active galactice nuclei jet plasmoids (diffusive case)](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/comparison/Comparison%20in%20compact%20astrophysical%20sources%201e17m%20trajectory.ipynb)
- [Ultra-high-energy cosmic rays in intergalactic magnetic fields](https://gitlab.ruhr-uni-bochum.de/reichp2y/proppy/-/blob/master/comparison/Comparison%20of%20UHECRs%20in%20IGMFs%204e27m.ipynb)

#### Basic example
```
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import proppy as pp
```
The software is modular. After initializing the simulation, we can add the indvidual modules as needed.

**Important modules are:**
- source
- magnetic field
- propagator
- observer

Each module contains special classes for different use cases. In the first tutorial, we only use basic classes to demonstrate the basic usage of the propagation software. In the following tutorials, each module with its special classes will be explained in detail. For example, there exist many different observer types such as time evolution observers (lin or log spcaing) and sphercial observers. 

```
sim = pp.Simulation()
```

#### Source
First, we have to define a source of particles that we want to propagate. The simples source type is the point source that emmits particles isotropically. The only user-specifyed parameters are:
- **energy**: Energy of the particles in eV.
- **source posistion**: The position of the point source.
- **number of particles**: Number of particles that should be emitted from this source.

The source can be easily added to the simulation. Afterwards, calling the description on the source of the simulation prints all relevant information and the values of the source parameters.

```
nr_particles = 10**3
source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
energy = 10**15 # eV

source = pp.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
sim.add_source(source)
sim.source.get_description()
```

#### Propagator
Propagates particles via a correlated rrandom walk. The overall behaviour is governed by a generalized telegraph equation.

Here, we use isotropic diffusion for simplicity that propagates particles using a correlated random walk in Cartesian coordinates. Isotropic diffusion is given when the turbulence is isotropic and there is no background field. In the other tutorials also the anisotropic diffusion is used. 

The only user-specifyed parameters for the simple case of isotropic diffusion are:
- **mean-free paths**: The mean-free paths λ​ can be derived from the diffusion coefficients κ​ via λ​=3κ​/c. Therfore, the diagonal elements of the diagonalized diffusion tensor are needed to determine the characteristics of the transport. For isotropic diffusion all diffusion coefficients are equal​. Typical free-mean paths of charged particles in plasmoids in AGN jets are 10^12m (see Reichherzer et al. (2021)).
- **number of steps**: The number of simulation steps for each individual particle.
- **step size**: Size of an individual step. Together with the parameter number of steps, the step size determines the trajectory length of the particles. 

The propagator can be easily added to the simulation. Afterwards, calling the description on the propagator of the simulation prints all relevant information and the values of the propagation parameters.

```
nr_steps = 10**4
step_size = 0.5*10**10 # [m]
diffusion_coefficient = 5*10**18 # [m^2/s]
speed_of_light = 3*10**8 # [m/s]
mfp_iso = 3*diffusion_coefficient/speed_of_light
mfp = np.array([10**11, 10**11, 10**11], dtype=np.float32)  # [m]
rms = 1 # Gaus
magnetic_field = pp.OrderedBackgroundField(rms, [0,0,1]).magnetic_field

propagator = pp.IsotropicPropagator( mfp, nr_steps, step_size)
sim.add_propagator(propagator)
sim.propagator.get_description()
```

#### Observer

The Observer determines during the simulation when and what data to write out (observe).In each simulation step, the current particle state is evaluated by the Observer to check if one of the observing contions is satisfyed. The conditions to observe can be based 
on the time (-> step) or the coordinates of the particle. The conditions to observe can be based on the time (or step) or the coordinates of the particle.

step number [unique_steps] -> time (TimeEvolutionObservers)radius of observer sphere [shperes] -> sphere around source (SphericalObservers)cartesian coordinates [box_dimensions] -> box around source (BoxObserver) (not yet implemented)
All special observer will create an Observer object and specify the relevant parameters
for the observation conditions (unique_steps, shperes, box_dimensions)

```
substeps = [False, False, True] # observe only steps (no substeps)
min_step = 1
max_step = nr_steps
nr_obs_steps = 200

observer = pp.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)

sim.add_observer(observer)
sim.observer.get_description()
```

#### Run simulation
```
sim.run_simulation()
sim.save_data('data/data_tut_1')
```

#### Analyze Statistics

As the particles propagate via a random walk, statistical properties of many particles are interesting, such as the diffusion coefficients and particle distributions. These quantities can be compared to analytical predictions. Running diffusion coefficients are computed with the mean-squared displacement method: κ(t) = <Δx²>/2t

```
df = pd.read_pickle("data/data_tut_1.pkl")
df_time_evolution_observer = df.loc[df['radius'] == -1.0]

sta = pp.Statistics(df_time_evolution_observer)
isotropic = True # diffusion is isotropic
errors = False # don't show error bars
df_kappas = sta.plot_diffusion_coefficients(isotropic, errors)
```
<img src="https://gitlab.ruhr-uni-bochum.de/tp4/diffusion/rwpropa/-/wikis/uploads/ba57a8da0e991d66291b0e6d64da11b5/kappa.png" alt="running diffusion coefficient" width="500"/>

# Scientific usage
In principle, this code is used wherever other propagation codes such as [CRPropa](https://crpropa.github.io/CRPropa3/), [DRAGON](https://github.com/cosmicrays/DRAGON), [GALPROP](https://galprop.stanford.edu/) are already in use. However, the advantages of RWPropa are especially in the good performance and the correct description also of the initial propagation, which is not possible for pure diffusive propagation approaches. Therefore, this code is predestined for the simulation of the particle transport in compact objects, where the initial, ballistic transport phase is often relevant (the importance depends on the scales and the energies).

As an example the propagation of charged particles in blobs of blazar jets is discussed in the following. Due to the high performance and the good statistical description even at early times, the software is excellently suited for the calculation of escape times of charged particles from certain zones (blob in the example), which in turn are required in (semi)analytical calculations. Also particle distributions and arrival times can be simulated efficiently. 

Contribution on [HEASA 2021 conference](https://indico.cern.ch/event/1037017/):
![poster](https://gitlab.ruhr-uni-bochum.de/tp4/diffusion/rwpropa/-/wikis/uploads/d7846a634cb0e85114f39f012663bcb4/poster_small.png)


# Documentation
For more detailed documentation of all parameters, return types, and behaviors of the above modules, please refer to the in-code documentation that is present in the header of each file, class and function in the software. In the following the most improtant modules are presented.

#### Simulation module
This module provides the Simulation class where all other modules are defined and set for the simulation, such as the source, the observer, and the propagator.

#### Source module
There are different sources available that can be customized by the user. 
The source specifies the initial state of the particles.

##### - PointSourceOriented
A point source that emits particles into a user-defined direction.

All particles start from a single point defined by the source position in the user-defined direction. All particles have the exact same state in the beginning.

##### - PointSourceIsotropicPhi
A point source that emits particles isotropically in phi.

All particles start from a single point defined by the source position in the user-defined direction. All particles have the exact same state in the beginning, except for the direction, which is isotropic in phi.

##### - PointSourceIsotropic
A point source that emits particles isotropically.

All particles start from a single point defined by the source position in the user-defined direction. All particles have the exact same state in the beginning, except for the direction, which is isotropic.

#### Observer module
The Observer determines during the simulation when and what data to write out (observe).

In each simulation step, the current particle state is evaluated by the Observer to check
if one of the observing conditions is satisfied. The conditions to observe can be based 
on the time (-> step) or the coordinates of the particle.


##### - ObserverAllSteps
Observes particles in all propagation steps.

##### - TimeEvolutionObserverLog
Observes particles at the user specified step numbers.
The user only gives the minimum, the maximum and the total step numbers. The TimeEvolutionObserverLog computes the list (logarithmically).

##### - TimeEvolutionObserverLin
Observes particles at the user specified step numbers.
The user only gives the minimum, the maximum and the total step numbers. The TimeEvolutionObserverLin computes the list (linearly).

##### - TimeEvolutionObserver
Observes particles at the user specified step numbers. The user passes the list of steps to the TimeEvolutionObserver.

#### Propagator module

_Background information:_
Current simulation codes propagate charged particles through magnetic fields either by ballistic codes solving the Lorentz force at each step or by using particle distributions arising from the solution of the classical diffusion equation. Whereas the former method provides an accurate description of the particle trajectories using small step sizes, this method is computationally expensive, especially for low particle energies. The method based on the diffusion equation, on the other hand, is extremely fast since the particle distribution is analytically given at each time point, but by the nature of the method, particles can only be described statistically. Even using quasi-particles, the individual trajectories are useless. For applications in which statistical statements and averaging over many particles are intended, this method is preferred in many areas of astronomy because of the short simulation times.

It is important to note, however, that the diffusion equation guarantees an adequate description of the particles only in the limit of infinitely large times. Numerically, however, this approach can also be used starting from times for which a diffusive behavior occurs. This typically happens after a simulation time which is in the order of magnitude of the mean free path length. 

Consequently, the use of propagation codes based on the diffusion equation is not suitable in compact objects. Examples include AGNs, supernovas, pulsars, etc.

_Implementation:_
In analogy to existing diffusion codes, the following routine requires the information of the diffusion coefficients. Here, however, we do not use the solution of the diffusion equation since this only provides a particle distribution and makes statistical statements about individual particle trajectories. Instead, we propagate particles according to the two-step propagation routine. 
- For **strong turbulence levels** we use the correlated random walk (CRW) in Cartesian coordinates, and
- For **weak turbulence levels** the CRW in cylindrical coordinates.

#### Magnetic field module
Ordered magnetic fields determine the directions of the parallel and perpendicular diffusion coefficients of the diffusion tensor. The parallel diffusion coefficient is used for the parallel direction to the ordered background magnetic field. 
The magnetic fields are used by the propagation class to determine the orientation of the field and how to align the diffusion tensor. 

##### - OrderedBackgroundField
The user can specify the root-mean-square (rms) value and the direction of the ordered magnetic field. The magnetic field is static and points everywhere in the specified direction with the rms value that is specified.

##### - DefaultBackgroundField
The user can only specify the root-mean-square (rms) value of the ordered magnetic field. 
The magnetic field is static and points everywhere along the z-axis with the rms value that is specified.


#### Trajectory module
Module to visualize simulated particle trajectories.

The simulated particles can be visualized here by showing their trajectories
or other properties.

After loading the simulation data, the simulated particles can be visualized here by showing their trajectories or other properties.

#### Statistics module
Scripts to compute statistical quantities of the particles.

As the particles propagate via a random walk, statistical properties of many particles are interesting, such as the diffusion coefficients and particle distributions. These quantities can be compared to analytical predictions.

After loading the simulation data, distributions of particles and running diffusion coefficients can be plotted.

The routine on how to compute diffusion coefficients is described in:

[1] Reichherzer, P., Becker Tjus, J., Zweibel, E. G., Merten, L., and Pueschel, M. J., 
“Turbulence-level dependence of cosmic ray parallel diffusion”, Monthly Notices of the Royal Astronomical Society, vol. 498, no. 4, pp. 5051–5064, 2020. doi:10.1093/mnras/staa2533.


# Community guidelines
### Contribute

Contributions are always welcome! Please fork the repository and open a new pull request for any new features.

After making changes, make sure everything works by running

`python3 -m unittest tests.py`

Please also add a new unittest for testing your new development if appropriate. 

### Support
Report issues or problems with the software or seek support by open new issues. 

doi: 10.5281/zenodo.5525925
