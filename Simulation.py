
import matplotlib.pyplot as plt
import random
import numpy as np
from numba import jit, int32, float32, types, typed
from numba.typed import List
from numba.experimental import jitclass

class Simulation():
    def __init__(self):
        self.particles
        self.propagation_mode
        self.time

    def addParticles(self, particles):
        self.particles = particles

    @jit(nopython=True)
    def runSimulation(self):
        kappa_perp = []
        kappa_para = []
        x = [0]
        y = [0]
        for i, t in enumerate(self.time):
            kappa_perp_sum = 0
            kappa_para_sum = 0
            for p in self.particles:
                p.move()
                kappa_perp_sum = kappa_perp_sum + (p.kappa(0)+p.kappa(1))/2
                kappa_para_sum = kappa_para_sum + p.kappa(2)
                if i == len(self.time)-1:
                    x.append(p.pos[0])
                    y.append(p.pos[1])
            
            kappa_para.append(kappa_para_sum/len(self.particles))
            kappa_perp.append(kappa_perp_sum/len(self.particles))

        return [kappa_para, kappa_perp, x, y]
