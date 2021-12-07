# run script via "python3 -m unittest tests.py"
import unittest
import numpy as np
import pandas as pd
import os
os.chdir('..')
import rwpropa as rw


class TestSource(unittest.TestCase):

    def test_point_source_oriented(self):
        print('\n----------------------------------')
        print('-> unit_test_point_source_oriented')
        energy = 10**10 #eV
        pos = [0,0,0]
        nr_particles = 10**2
        pitch_angle = 2*np.pi * 54.74/360
        phi = np.pi/4.0
        point_source_oriented = rw.PointSourceOriented(energy, pos, nr_particles, pitch_angle, phi)
        # test if the number of generated particles is correct
        self.assertEqual(nr_particles, len(point_source_oriented.particles))
        # test if the energy of a generated particle is correct
        self.assertEqual(energy, point_source_oriented.particles[0].ps.energy)
        # check the position of the particles
        expected_pos = list(np.array(pos, dtype=np.float32))
        actual_pos = point_source_oriented.particles[0].ps.pos.tolist()
        self.assertListEqual(expected_pos, actual_pos)
        # check the direction of the particles
        precision = 4
        expected_direction = list(np.around(np.array([
            np.cos(phi)*np.sin(pitch_angle), 
            np.sin(phi)*np.sin(pitch_angle), 
            np.cos(pitch_angle)
        ], dtype=np.float32), precision))
        actual_direction = list(np.around(np.array(point_source_oriented.particles[0].ps.dir), precision))
        self.assertListEqual(expected_direction, actual_direction)

    
    def test_point_source_isotropic(self):
        print('\n----------------------------------')
        print('-> unit_test_point_source_isotropic')
        energy = 10**10 #eV
        pos = [0,0,0]
        nr_particles = 10**2
        point_source_oriented = rw.PointSourceIsotropic(energy, pos, nr_particles)
        # test if the number of generated particles is correct
        self.assertEqual(nr_particles, len(point_source_oriented.particles))
        # test if the energy of a generated particle is correct
        self.assertEqual(energy, point_source_oriented.particles[0].ps.energy)
        # check the position of the particles
        expected_pos = list(np.array(pos, dtype=np.float32))
        actual_pos = point_source_oriented.particles[0].ps.pos.tolist()
        self.assertListEqual(expected_pos, actual_pos)



class TestObserver(unittest.TestCase):

    def test_time_evolution_observer_unit(self):
        print('\n----------------------------------')
        print('-> unit_test_time_evolution_observer')
        sim = rw.Simulation()
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = 100
        nr_obs_steps = 10
        observer = rw.TimeEvolutionObserverLin(min_step, max_step, nr_obs_steps, substeps)
        sim.add_observer(observer)
        self.assertEqual(min_step, sim.observer.steps[0])
        self.assertEqual(max_step, sim.observer.steps[-1])



class TestIntegration(unittest.TestCase):

    def test_basic_propagation_isotropic_source(self):
        print('\n----------------------------------')
        print('-> integration_test_basic_propagation_isotropic_source')

        sim = rw.Simulation()

        # adding a particle source
        energy = 10**10 #eV
        pos = [0,0,0]
        nr_particles = 10**1
        source = rw.PointSourceIsotropic(energy, pos, nr_particles)
        sim.add_source(source)
        
        # adding a propagator to simulation
        nr_steps = 10**3
        step_size = 0.5*10**10 # [m]
        mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
        rms = 1 # Gaus
        magnetic_field = rw.OrderedBackgroundField(rms, [0,0,1]).magnetic_field
        propagator = rw.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
        sim.add_propagator(propagator)

        # adding a TimeEvolutionObserver
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = nr_steps
        nr_obs_steps = 30
        observer = rw.TimeEvolutionObserverLin(min_step, max_step, nr_obs_steps, substeps)
        sim.add_observer(observer)

        # simulate
        sim.run_simulation()
        df = pd.DataFrame(sim.data[1:])

        # check if the number of observations is correct
        self.assertEqual(nr_obs_steps*nr_particles, len(df[0]))

        # check if the step_number of the first observed step is correct
        self.assertEqual(min_step, df[1].tolist()[0])

        # check if the step_number of the last observed step is correct
        self.assertEqual(max_step, df[1].tolist()[-1])


    def test_basic_propagation_oriented_source(self):
        print('\n----------------------------------')
        print('-> integration_test_basic_propagation_oriented_source')

        sim = rw.Simulation()

        # adding a particle source
        energy = 10**10 #eV
        pos = [0,0,0]
        nr_particles = 10**1
        pitch_angle = 2*np.pi * 54.74/360
        phi = np.pi/4.0
        source = rw.PointSourceOriented(energy, pos, nr_particles, pitch_angle, phi)
        sim.add_source(source)

        # adding a propagator to simulation
        nr_steps = 10**3
        step_size = 0.5*10**10 # [m]
        mfp = np.array([2.13*10**12/2.0, 2.13*10**12/2.0, 2.1078*10**12], dtype=np.float32)  # [m]
        rms = 1 # Gaus
        magnetic_field = rw.OrderedBackgroundField(rms, [0,0,1]).magnetic_field
        propagator = rw.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
        sim.add_propagator(propagator)

        # adding a TimeEvolutionObserver
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = nr_steps
        nr_obs_steps = 30
        observer = rw.TimeEvolutionObserverLin(min_step, max_step, nr_obs_steps, substeps)
        sim.add_observer(observer)

        # simulate
        sim.run_simulation()
        df = pd.DataFrame(sim.data[1:])

        # check if the number of observations is correct
        self.assertEqual(nr_obs_steps*nr_particles, len(df[0]))

        # check if the step_number of the first observed step is correct
        self.assertEqual(min_step, df[1].tolist()[0])

        # check if the step_number of the last observed step is correct
        self.assertEqual(max_step, df[1].tolist()[-1])

    
    def test_isotropic_diffusion_coefficient(self):
        # Diffusion coefficients are a statistical description of the transport properties of charged particles in turbulent magnetic fields. Monte Carlo simulations of charged particles can be used to calculate diffusion coefficients. In this integration test, we determine the diffusion coefficient for a parameter configuration for which we know the diffusion coefficients. For this, we need to build our simulation as follows. After initializing the simulation, we need add the indvidual modules:
        # - source
        # - magnetic field
        # - propagator
        # - observer
        print('\n----------------------------------')
        print('-> test_isotropic_diffusion_coefficient')

        sim = rw.Simulation()

        # adding a particle source
        nr_particles = 10**3
        source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        energy = 10**15 # eV
        source = rw.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
        sim.add_source(source)

        # adding a propagator to simulation
        nr_steps = 10**3
        step_size = 0.2*10**10 # [m]
        diffusion_coefficient = 5*10**18 # [m^2/s]
        speed_of_light = 3*10**8 # [m/s]
        mfp_iso = 3*diffusion_coefficient/speed_of_light
        mfp = np.array([mfp_iso, mfp_iso, mfp_iso], dtype=np.float32)  # [m]
        propagator = rw.IsotropicPropagator(mfp, nr_steps, step_size)
        sim.add_propagator(propagator)

        # adding a TimeEvolutionObserver
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = nr_steps
        nr_obs_steps = 600
        observer = rw.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)
        sim.add_observer(observer)

        # simulate
        print('Simulation can take some time ...')
        sim.run_simulation()
        df = pd.DataFrame(sim.data[1:])

        # check if the step_number of the last observed step is correct
        self.assertEqual(max_step, df[1].tolist()[-1])

        # get diffusion coefficients by using the statistics class
        df.columns = observer.get_column_names()
        sta = rw.Statistics(df)
        isotropic = True # diffusion is isotropic
        errors = False # don't show error bars
        df_kappas = sta.plot_diffusion_coefficients(isotropic, errors, None)
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
        self.assertTrue(diffusion_coefficient*0.8 <= kappa_xx <= diffusion_coefficient*1.2)
        # test if kappa_yy is in expected range
        self.assertTrue(diffusion_coefficient*0.8 <= kappa_yy <= diffusion_coefficient*1.2)
        # test if kappa_zz is in expected range
        self.assertTrue(diffusion_coefficient*0.8 <= kappa_zz <= diffusion_coefficient*1.2)
        # test if kappa is in expected range
        kappa = np.mean(np.array([kappa_xx, kappa_yy, kappa_zz]))
        kappa_err = np.std(np.array([kappa_xx, kappa_yy, kappa_zz]))
        print('kappa:', f"{kappa:.3}", 'm²/s', '+-', f"{kappa_err:.3}", 'm²/s')
        self.assertTrue(diffusion_coefficient*0.9 <= kappa <= diffusion_coefficient*1.1)


    def test_anisotropic_diffusion_coefficient(self):
        # Diffusion coefficients are a statistical description of the transport properties of charged particles in turbulent magnetic fields. Monte Carlo simulations of charged particles can be used to calculate diffusion coefficients. In this integration test, we determine the anisotropic diffusion tensor for a parameter configuration for which we know the diffusion coefficients. For this, we need to build our simulation as follows. After initializing the simulation, we need add the indvidual modules:
        # - source
        # - magnetic field
        # - propagator
        # - observer
        print('\n----------------------------------')
        print('-> test_anisotropic_diffusion_coefficient')

        sim = rw.Simulation()

        # adding a particle source
        nr_particles = 3*10**2
        source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        energy = 3*10**15 # eV

        source = rw.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
        sim.add_source(source)

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
        magnetic_field = rw.OrderedBackgroundField(rms, [0,0,1]).magnetic_field

        propagator = rw.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
        sim.add_propagator(propagator)

        # adding a TimeEvolutionObserver
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = nr_steps
        nr_obs_steps = 600

        observer = rw.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)

        sim.add_observer(observer)

        # simulate
        print('Simulation can take some time ...')
        sim.run_simulation()
        df = pd.DataFrame(sim.data[1:])

        # check if the step_number of the last observed step is correct
        self.assertEqual(max_step, df[1].tolist()[-1])

        # get diffusion coefficients by using the statistics class
        df.columns = observer.get_column_names()
        sta = rw.Statistics(df)
        isotropic = True # diffusion is isotropic
        errors = False # don't show error bars
        df_kappas = sta.plot_diffusion_coefficients(isotropic, errors, None)
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
        self.assertTrue(diffusion_coefficient_perp*(1.-uncertainty) <= kappa_xx <= diffusion_coefficient_perp*(1.+uncertainty))
        # test if kappa_yy is in expected range
        self.assertTrue(diffusion_coefficient_perp*(1.-uncertainty) <= kappa_yy <= diffusion_coefficient_perp*(1.+uncertainty))
        # test if kappa_zz is in expected range
        self.assertTrue(diffusion_coefficient_para*(1.-uncertainty) <= kappa_zz <= diffusion_coefficient_para*(1.+uncertainty))
        # test if kappa is in expected range
        kappa_perp = np.mean(np.array([kappa_xx, kappa_yy]))
        kappa_perp_err = np.std(np.array([kappa_xx, kappa_yy]))
        print('kappa:', f"{kappa_perp:.3}", 'm²/s', '+-', f"{kappa_perp_err:.3}", 'm²/s')
        self.assertTrue(diffusion_coefficient_perp*(1.-uncertainty) <= kappa_perp <= diffusion_coefficient_perp*(1.+uncertainty))


    def test_simple_anisotropic_diffusion_coefficient(self):
        # Similar to test_anisotropic_diffusion_coefficient, but much faster as only tested if the parallel 
        # diffusion coefficient is larger than the perpendicular diffusion coefficient.
        print('\n----------------------------------')
        print('-> test_simple_anisotropic_diffusion_coefficient')

        sim = rw.Simulation()

        # adding a particle source
        nr_particles = 10**2
        source_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        energy = 3*10**15 # eV

        source = rw.PointSourceIsotropicPhi(energy, source_pos, nr_particles)
        sim.add_source(source)

        # adding a propagator to simulation
        nr_steps = 1*10**4
        step_size = 0.5*10**10 # [m]
        speed_of_light = 3*10**8 # [m/s]
        diffusion_coefficient_perp = 1*10**17 # [m^2/s]
        diffusion_coefficient_para = 1*10**20 # [m^2/s]
        mfp_perp = 3*diffusion_coefficient_perp/speed_of_light*2
        mfp_para = 3*diffusion_coefficient_para/speed_of_light
        mfp = np.array([mfp_perp, mfp_perp, mfp_para], dtype=np.float32)
        rms = 1 # Gaus
        magnetic_field = rw.OrderedBackgroundField(rms, [0,0,1]).magnetic_field

        propagator = rw.AnisotropicPropagator(magnetic_field, mfp, nr_steps, step_size)
        sim.add_propagator(propagator)

        # adding a TimeEvolutionObserver
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = nr_steps
        nr_obs_steps = 200

        observer = rw.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)

        sim.add_observer(observer)

        # simulate
        print('Simulation can take some time ...')
        sim.run_simulation()
        df = pd.DataFrame(sim.data[1:])

        # check if the step_number of the last observed step is correct
        self.assertEqual(max_step, df[1].tolist()[-1])

        # get diffusion coefficients by using the statistics class
        df.columns = observer.get_column_names()
        sta = rw.Statistics(df)
        isotropic = True # diffusion is isotropic
        errors = False # don't show error bars
        df_kappas = sta.plot_diffusion_coefficients(isotropic, errors, None)
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
        uncertainty = 0.6
        # test if kappa_xx is in expected range
        self.assertTrue(diffusion_coefficient_perp*(1.-uncertainty) <= kappa_xx <= diffusion_coefficient_perp*(1.+uncertainty))
        # test if kappa_yy is in expected range
        self.assertTrue(diffusion_coefficient_perp*(1.-uncertainty) <= kappa_yy <= diffusion_coefficient_perp*(1.+uncertainty))
        # test if kappa_zz is in expected range
        self.assertTrue(diffusion_coefficient_para*(1.-uncertainty) <= kappa_zz <= diffusion_coefficient_para*(1.+uncertainty))
        # test if kappa is in expected range
        kappa_perp = np.mean(np.array([kappa_xx, kappa_yy]))
        kappa_perp_err = np.std(np.array([kappa_xx, kappa_yy]))
        print('kappa:', f"{kappa_perp:.3}", 'm²/s', '+-', f"{kappa_perp_err:.3}", 'm²/s')
        self.assertTrue(diffusion_coefficient_perp*(1.-uncertainty) <= kappa_perp <= diffusion_coefficient_perp*(1.+uncertainty))