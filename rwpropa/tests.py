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
        nr_steps = 10**4
        step_size = 0.2*10**10 # [m]
        diffusion_coefficient = 5*10**18 # [m^2/s]
        speed_of_light = 3*10**8 # [m/s]
        mfp_iso = 3*diffusion_coefficient/speed_of_light
        mfp = np.array([mfp_iso, mfp_iso, mfp_iso], dtype=np.float32)  # [m]

        # adding a TimeEvolutionObserver
        propagator = rw.IsotropicPropagator(mfp, nr_steps, step_size)
        sim.add_propagator(propagator)
        substeps = [False, False, True] # observe only steps (no substeps)
        min_step = 1
        max_step = nr_steps
        nr_obs_steps = 600
        observer = rw.TimeEvolutionObserverLog(min_step, max_step, nr_obs_steps, substeps)
        sim.add_observer(observer)

        # simulate
        sim.run_simulation()
        df = pd.DataFrame(sim.data[1:])

        # check if the step_number of the last observed step is correct
        self.assertEqual(max_step, df[1].tolist()[-1])