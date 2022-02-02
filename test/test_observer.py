import numpy as np
import pandas as pd
import proppy.simulation as simulation
import proppy.observer as observer
import proppy.source as source
import proppy.propagator as propagator


def test_time_evolution_observer_unit():
        print('\n----------------------------------')
        print('-> unit_test_time_evolution_observer')
        sim = simulation.Simulation()
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = 100
        nr_obs_steps = 10
        obs = observer.TimeEvolutionObserverLin(min_step, max_step, nr_obs_steps, substeps)
        sim.add_observer(obs)
        assert min_step == sim.observer.steps[0]
        assert max_step == sim.observer.steps[-1]


def test_spherical_observer_unit():
    print('\n----------------------------------')
    print('-> unit_test_spherical_observer')
    sim = simulation.Simulation()
    nr_particles = 1
    source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    energy = 10**12 # eV
    src = source.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
    sim.add_source(src)

    nr_steps = 1*10**5
    step_size = 1.0*10**9 # [m]
    diffusion_coefficient = 1.5*10**20 # [m^2/s]
    speed_of_light = 3*10**8 # [m/s]
    mfp_iso = 3*diffusion_coefficient/speed_of_light
    mfp = np.array([mfp_iso, mfp_iso, mfp_iso], dtype=np.float32)  # [m]
    prop = propagator.IsotropicPropagator(mfp, nr_steps, step_size)
    sim.add_propagator(prop)

    substeps = [False, False, True] # observe only steps (no substeps)
    sphere = 10**10 # [m]
    spheres = [sphere]
    obs = observer.SphericalObserver(substeps, spheres, on_detection_deactivate=True)
    sim.add_observer(obs)

    sim.run_simulation()
    df = pd.DataFrame(sim.data[1:])

    # check if the number of observations is correct -> only observe all particles once
    assert nr_particles == len(df[0])

    # check if the trajcetory lenth of the observed paarticle is at least the sphere radius
    trajectory_length = df[2].tolist()[0]
    assert sphere <= trajectory_length
    # check if the trajcetory lenth of the observed paarticle is not much larger than the sphere radius
    # because we assume ballistic propagation phase until its observation given the large diffusion coefficient
    assert trajectory_length <= 2*sphere 