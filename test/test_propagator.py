import numpy as np
import pandas as pd
import proppy.simulation as simulation
import proppy.observer as observer
import proppy.source as source
import proppy.propagator as propagator
import proppy.magnetic_field as mag_field


def test_basic_propagation_isotropic_source():
    print('\n----------------------------------')
    print('-> integration_test_basic_propagation_isotropic_source')

    sim = simulation.Simulation()

    # adding a particle source
    energy = 10**10 #eV
    pos = [0,0,0]
    nr_particles = 10**1
    src = source.PointSourceIsotropic(energy, pos, nr_particles)
    sim.add_source(src)
    
    # adding a propagator to simulation
    nr_steps = 10**3
    step_size = 0.5*10**10 # [m]
    mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
    rms = 1 # Gaus
    magnetic_field = mag_field.OrderedBackgroundField(rms, [0,0,1]).magnetic_field
    prop = propagator.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim.add_propagator(prop)

    # adding a TimeEvolutionObserver
    substeps = [False, False, True] # observe only steps (no substeps)
    min_step = 1
    max_step = nr_steps
    nr_obs_steps = 30
    obs = observer.TimeEvolutionObserverLin(min_step, max_step, nr_obs_steps, substeps)
    sim.add_observer(obs)

    # simulate
    sim.run_simulation()
    df = pd.DataFrame(sim.data[1:])

    # check if the number of observations is correct
    assert nr_obs_steps*nr_particles == len(df[0])

    # check if the step_number of the first observed step is correct
    assert min_step == df[1].tolist()[0]

    # check if the step_number of the last observed step is correct
    assert max_step == df[1].tolist()[-1]


def test_basic_propagation_oriented_source():
    print('\n----------------------------------')
    print('-> integration_test_basic_propagation_oriented_source')

    sim = simulation.Simulation()

    # adding a particle source
    energy = 10**10 #eV
    pos = [0,0,0]
    nr_particles = 10**1
    pitch_angle = 2*np.pi * 54.74/360
    phi = np.pi/4.0
    src = source.PointSourceOriented(energy, pos, nr_particles, pitch_angle, phi)
    sim.add_source(src)

    # adding a propagator to simulation
    nr_steps = 10**3
    step_size = 0.5*10**10 # [m]
    mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
    rms = 1 # Gaus
    magnetic_field = mag_field.OrderedBackgroundField(rms, [0,0,1]).magnetic_field
    prop = propagator.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim.add_propagator(prop)

    # adding a TimeEvolutionObserver
    substeps = [False, False, True] # observe only steps (no substeps)
    min_step = 1
    max_step = nr_steps
    nr_obs_steps = 30
    obs = observer.TimeEvolutionObserverLin(min_step, max_step, nr_obs_steps, substeps)
    sim.add_observer(obs)

    # simulate
    sim.run_simulation()
    df = pd.DataFrame(sim.data[1:])

    # check if the number of observations is correct
    assert nr_obs_steps*nr_particles == len(df[0])

    # check if the step_number of the first observed step is correct
    assert min_step == df[1].tolist()[0]

    # check if the step_number of the last observed step is correct
    assert max_step == df[1].tolist()[-1]