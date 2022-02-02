import numpy as np
import pandas as pd
import proppy.simulation as simulation
import proppy.observer as observer
import proppy.source as source
import proppy.propagator as propagator
import proppy.magnetic_field as mag_field
import proppy.statistics as statistics

def test_isotropic_diffusion_coefficient():
    # Diffusion coefficients are a statistical description of the transport properties of charged particles in turbulent magnetic fields. Monte Carlo simulations of charged particles can be used to calculate diffusion coefficients. In this integration test, we determine the diffusion coefficient for a parameter configuration for which we know the diffusion coefficients. For this, we need to build our simulation as follows. After initializing the simulation, we need add the indvidual modules:
    # - source
    # - magnetic field
    # - propagator
    # - observer
    print('\n----------------------------------')
    print('-> test_isotropic_diffusion_coefficient')

    sim = simulation.Simulation()

    # adding a particle source
    nr_particles = 10**3
    source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    energy = 10**15 # eV
    src = source.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
    sim.add_source(src)

    # adding a propagator to simulation
    nr_steps = 10**3
    step_size = 0.2*10**10 # [m]
    diffusion_coefficient = 5*10**18 # [m^2/s]
    speed_of_light = 3*10**8 # [m/s]
    mfp_iso = 3*diffusion_coefficient/speed_of_light
    mfp = np.array([mfp_iso, mfp_iso, mfp_iso], dtype=np.float32)  # [m]
    prop = propagator.IsotropicPropagator(mfp, nr_steps, step_size)
    sim.add_propagator(prop)

    # adding a TimeEvolutionObserver
    substeps = [False, False, True] # observe only steps (no substeps)
    min_step = 1
    max_step = nr_steps
    nr_obs_steps = 600
    obs = observer.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)
    sim.add_observer(obs)

    # simulate
    print('Simulation can take some time ...')
    sim.run_simulation()
    df = pd.DataFrame(sim.data[1:])

    # check if the step_number of the last observed step is correct
    assert max_step == df[1].tolist()[-1]

    # get diffusion coefficients by using the statistics class
    df.columns = obs.get_column_names()
    sta = statistics.Statistics(df)
    isotropic = True # diffusion is isotropic
    errors = False # don't show error bars
    df_kappas = sta.plot_diffusion_coefficients(isotropic, errors, None, plot=False)
    print('input kappa:', f"{float(diffusion_coefficient):.3}", 'm²/s')
    n = 20
    kappa_xx = np.mean(df_kappas['kappa_xx'][-n:])
    kappa_yy = np.mean(df_kappas['kappa_yy'][-n:])
    kappa_zz = np.mean(df_kappas['kappa_zz'][-n:])
    kappa_xx_err = np.std(df_kappas['kappa_xx'][-n:])
    kappa_yy_err = np.std(df_kappas['kappa_yy'][-n:])
    kappa_zz_err = np.std(df_kappas['kappa_zz'][-n:])
    print('kappa_{xx}:', f"{kappa_xx:.3}", 'm²/s', '+-', f"{kappa_xx_err:.3}", 'm²/s')
    print('kappa_{yy}:', f"{kappa_yy:.3}", 'm²/s', '+-', f"{kappa_yy_err:.3}", 'm²/s')
    print('kappa_{zz}:', f"{kappa_zz:.3}", 'm²/s', '+-', f"{kappa_zz_err:.3}", 'm²/s')
    print('Note that there is an additional systematic error that can lead to minor deviations between theory and simulations given the limited particle trajectory length. When increasing the trajectory length the agreement improves, but the simulations take longer.')

    # test if kappa_xx is in expected range
    assert diffusion_coefficient*0.8 <= kappa_xx <= diffusion_coefficient*1.2
    # test if kappa_yy is in expected range
    assert diffusion_coefficient*0.8 <= kappa_yy <= diffusion_coefficient*1.2
    # test if kappa_zz is in expected range
    assert diffusion_coefficient*0.8 <= kappa_zz <= diffusion_coefficient*1.2
    # test if kappa is in expected range
    kappa = np.mean(np.array([kappa_xx, kappa_yy, kappa_zz]))
    kappa_err = np.std(np.array([kappa_xx, kappa_yy, kappa_zz]))
    print('kappa:', f"{kappa:.3}", 'm²/s', '+-', f"{kappa_err:.3}", 'm²/s')
    assert diffusion_coefficient*0.9 <= kappa <= diffusion_coefficient*1.1


def test_anisotropic_diffusion_coefficient():
    # Diffusion coefficients are a statistical description of the transport properties of charged particles in turbulent magnetic fields. Monte Carlo simulations of charged particles can be used to calculate diffusion coefficients. In this integration test, we determine the anisotropic diffusion tensor for a parameter configuration for which we know the diffusion coefficients. For this, we need to build our simulation as follows. After initializing the simulation, we need add the indvidual modules:
    # - source
    # - magnetic field
    # - propagator
    # - observer
    print('\n----------------------------------')
    print('-> test_anisotropic_diffusion_coefficient')

    sim = simulation.Simulation()

    # adding a particle source
    nr_particles = 3*10**2
    source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    energy = 3*10**15 # eV

    src = source.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
    sim.add_source(src)

    # adding a propagator to simulation
    nr_steps = 1*10**5
    step_size = 0.2*10**10 # [m]
    speed_of_light = 3*10**8 # [m/s]
    diffusion_coefficient_perp = 1.3*10**18 # [m^2/s]
    diffusion_coefficient_para = 1.4*10**20 # [m^2/s]
    mfp_perp = 3*diffusion_coefficient_perp/speed_of_light*2
    mfp_para = 3*diffusion_coefficient_para/speed_of_light
    mfp = np.array([mfp_perp, mfp_perp, mfp_para], dtype=np.float32)
    rms = 1 # Gaus
    magnetic_field = mag_field.OrderedBackgroundField(rms, [0,0,1]).magnetic_field

    prop = propagator.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim.add_propagator(prop)

    # adding a TimeEvolutionObserver
    substeps = [False, False, True] # observe only steps (no substeps)
    min_step = 1
    max_step = nr_steps
    nr_obs_steps = 600

    obs = observer.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)
    sim.add_observer(obs)

    # simulate
    print('Simulation can take some time ...')
    sim.run_simulation()
    df = pd.DataFrame(sim.data[1:])

    # check if the step_number of the last observed step is correct
    assert max_step == df[1].tolist()[-1]

    # get diffusion coefficients by using the statistics class
    df.columns = obs.get_column_names()
    sta = statistics.Statistics(df)
    isotropic = True # diffusion is isotropic
    errors = False # don't show error bars
    df_kappas = sta.plot_diffusion_coefficients(isotropic, errors, None, plot=False)
    print('input kappa_perp:', f"{float(diffusion_coefficient_perp):.3}", 'm²/s')
    print('input kappa_para:', f"{float(diffusion_coefficient_para):.3}", 'm²/s')
    n = 20
    kappa_xx = np.mean(df_kappas['kappa_xx'][-n:])
    kappa_yy = np.mean(df_kappas['kappa_yy'][-n:])
    kappa_zz = np.mean(df_kappas['kappa_zz'][-n:])
    kappa_xx_err = np.std(df_kappas['kappa_xx'][-n:])
    kappa_yy_err = np.std(df_kappas['kappa_yy'][-n:])
    kappa_zz_err = np.std(df_kappas['kappa_zz'][-n:])
    print('kappa_{xx}:', f"{kappa_xx:.3}", 'm²/s', '+-', f"{kappa_xx_err:.3}", 'm²/s')
    print('kappa_{yy}:', f"{kappa_yy:.3}", 'm²/s', '+-', f"{kappa_yy_err:.3}", 'm²/s')
    print('kappa_{zz}:', f"{kappa_zz:.3}", 'm²/s', '+-', f"{kappa_zz_err:.3}", 'm²/s')
    print('Note that there is an additional systematic error that can lead to minor deviations between theory and simulations given the limited particle trajectory length. When increasing the trajectory length the agreement improves, but the simulations take longer.')

    # Given that we decrease the statistics, we expect some statistical fluctuation between theory and anayltics.
    # In tutorial 2, this example is shown with better statistics and good agreement.
    uncertainty = 0.3

    # test if kappa_xx is in expected range
    assert diffusion_coefficient_perp*(1.-uncertainty) <= kappa_xx <= diffusion_coefficient_perp*(1.+uncertainty)

    # test if kappa_yy is in expected range
    assert diffusion_coefficient_perp*(1.-uncertainty) <= kappa_yy <= diffusion_coefficient_perp*(1.+uncertainty)

    # test if kappa_zz is in expected range
    assert diffusion_coefficient_para*(1.-uncertainty) <= kappa_zz <= diffusion_coefficient_para*(1.+uncertainty)
    # test if kappa is in expected range
    kappa_perp = np.mean(np.array([kappa_xx, kappa_yy]))
    kappa_perp_err = np.std(np.array([kappa_xx, kappa_yy]))
    print('kappa:', f"{kappa_perp:.3}", 'm²/s', '+-', f"{kappa_perp_err:.3}", 'm²/s')
    assert diffusion_coefficient_perp*(1.-uncertainty) <= kappa_perp <= diffusion_coefficient_perp*(1.+uncertainty)


def test_simple_anisotropic_diffusion_coefficient():
    # Similar to test_anisotropic_diffusion_coefficient, but much faster as only tested if the parallel 
    # diffusion coefficient is larger than the perpendicular diffusion coefficient.
    print('\n----------------------------------')
    print('-> test_simple_anisotropic_diffusion_coefficient')

    sim = simulation.Simulation()

    # adding a particle source
    nr_particles = 10**2
    source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    energy = 3*10**15 # eV

    src = source.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
    sim.add_source(src)

    # adding a propagator to simulation
    nr_steps = 1*10**4
    step_size = 0.2*10**10 # [m]
    speed_of_light = 3*10**8 # [m/s]
    diffusion_coefficient_perp = 1*10**18 # [m^2/s]
    diffusion_coefficient_para = 1*10**20 # [m^2/s]
    mfp_perp = 3*diffusion_coefficient_perp/speed_of_light*2
    mfp_para = 3*diffusion_coefficient_para/speed_of_light
    mfp = np.array([mfp_perp, mfp_perp, mfp_para], dtype=np.float32)
    rms = 1 # Gaus
    magnetic_field = mag_field.OrderedBackgroundField(rms, [0,0,1]).magnetic_field

    prop = propagator.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim.add_propagator(prop)

    # adding a TimeEvolutionObserver
    substeps = [False, False, True] # observe only steps (no substeps)
    min_step = 1
    max_step = nr_steps
    nr_obs_steps = 200

    obs = observer.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)
    sim.add_observer(obs)

    # simulate
    print('Simulation can take some time ...')
    sim.run_simulation()
    df = pd.DataFrame(sim.data[1:])

    # check if the step_number of the last observed step is correct
    assert max_step == df[1].tolist()[-1]

    # get diffusion coefficients by using the statistics class
    df.columns = obs.get_column_names()
    sta = statistics.Statistics(df)
    isotropic = True # diffusion is isotropic
    errors = False # don't show error bars
    df_kappas = sta.plot_diffusion_coefficients(isotropic, errors, None, plot=False)
    print('input kappa_perp:', f"{float(diffusion_coefficient_perp):.3}", 'm²/s')
    print('input kappa_para:', f"{float(diffusion_coefficient_para):.3}", 'm²/s')
    n = 20
    kappa_xx = np.mean(df_kappas['kappa_xx'][-n:])
    kappa_yy = np.mean(df_kappas['kappa_yy'][-n:])
    kappa_zz = np.mean(df_kappas['kappa_zz'][-n:])
    kappa_xx_err = np.std(df_kappas['kappa_xx'][-n:])
    kappa_yy_err = np.std(df_kappas['kappa_yy'][-n:])
    kappa_zz_err = np.std(df_kappas['kappa_zz'][-n:])
    print('kappa_{xx}:', f"{kappa_xx:.3}", 'm²/s', '+-', f"{kappa_xx_err:.3}", 'm²/s')
    print('kappa_{yy}:', f"{kappa_yy:.3}", 'm²/s', '+-', f"{kappa_yy_err:.3}", 'm²/s')
    print('kappa_{zz}:', f"{kappa_zz:.3}", 'm²/s', '+-', f"{kappa_zz_err:.3}", 'm²/s')
    print('Note that there is an additional systematic error that can lead to minor deviations between theory and simulations given the limited particle trajectory length. When increasing the trajectory length the agreement improves, but the simulations take longer.')

    # Given that we decrease the statistics, we expect some statistical fluctuation between theory and anayltics.
    # test if kappa is in expected range
    kappa_perp = np.mean(np.array([kappa_xx, kappa_yy]))
    kappa_perp_err = np.std(np.array([kappa_xx, kappa_yy]))

    assert kappa_perp <= kappa_zz
    print('kappa:', f"{kappa_perp:.3}", 'm²/s', '+-', f"{kappa_perp_err:.3}", 'm²/s')
    


def test_ballistic_anisotropic_diffusion_coefficient():
    # See Reichherzer et al. (2021) for the difference between initial ballistic and final diffusive
    # propagation.
    print('\n----------------------------------')
    print('-> test_ballistic_anisotropic_diffusion_coefficient')

    sim = simulation.Simulation()

    # adding a particle source
    nr_particles = 1*10**2
    source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    energy = 3*10**15 # eV

    src = source.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
    sim.add_source(src)

    # adding a propagator to simulation
    nr_steps = 1*10**3
    step_size = 1*10**6 # [m]
    speed_of_light = 3*10**8 # [m/s]
    diffusion_coefficient_perp = 1.3*10**18 # [m^2/s]
    diffusion_coefficient_para = 1.4*10**20 # [m^2/s]
    mfp_perp = 3*diffusion_coefficient_perp/speed_of_light*2
    mfp_para = 3*diffusion_coefficient_para/speed_of_light
    mfp = np.array([mfp_perp, mfp_perp, mfp_para], dtype=np.float32)
    rms = 1 # Gaus
    magnetic_field = mag_field.OrderedBackgroundField(rms, [0,0,1]).magnetic_field

    prop = propagator.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
    sim.add_propagator(prop)

    # adding a TimeEvolutionObserver
    substeps = [False, False, True] # observe only steps (no substeps)
    min_step = 1
    max_step = nr_steps
    nr_obs_steps = 200

    obs = observer.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)
    sim.add_observer(obs)

    # simulate
    print('Simulation can take some time ...')
    sim.run_simulation()
    df = pd.DataFrame(sim.data[1:])

    # check if the step_number of the last observed step is correct
    assert max_step == df[1].tolist()[-1]

    # get diffusion coefficients by using the statistics class
    df.columns = obs.get_column_names()
    sta = statistics.Statistics(df)
    isotropic = True # diffusion is isotropic
    errors = False # don't show error bars
    df_kappas = sta.plot_diffusion_coefficients(isotropic, errors, None, plot=False)
    n = 20
    kappa_xx = np.mean(df_kappas['kappa_xx'][-n:])
    kappa_yy = np.mean(df_kappas['kappa_yy'][-n:])
    kappa_zz = np.mean(df_kappas['kappa_zz'][-n:])
    kappa_xx_err = np.std(df_kappas['kappa_xx'][-n:])
    kappa_yy_err = np.std(df_kappas['kappa_yy'][-n:])
    kappa_zz_err = np.std(df_kappas['kappa_zz'][-n:])
    print('kappa_{xx}:', f"{kappa_xx:.3}", 'm²/s', '+-', f"{kappa_xx_err:.3}", 'm²/s')
    print('kappa_{yy}:', f"{kappa_yy:.3}", 'm²/s', '+-', f"{kappa_yy_err:.3}", 'm²/s')
    print('kappa_{zz}:', f"{kappa_zz:.3}", 'm²/s', '+-', f"{kappa_zz_err:.3}", 'm²/s')
    print('Note that we are in the ballistic propagation phase.')

    kappa_perp = np.mean(np.array([kappa_xx, kappa_yy]))
    kappa_perp_err = np.std(np.array([kappa_xx, kappa_yy]))
    print('kappa:', f"{kappa_perp:.3}", 'm²/s', '+-', f"{kappa_perp_err:.3}", 'm²/s')
    # both diffusion coefficients should be similar in the ballistic phase:
    assert kappa_zz*0.9 <= kappa_perp <= kappa_zz*1.1