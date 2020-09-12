# -*- coding: utf-8 -*-
"""
Time: 2020/Aug/20

@author: Benchi Zhao
Python version: 3.7
Functionality: This module calculate the drag force (air resistance) when droplet moves in the air
Dependence: Need import 'input', all parameters of the experiment is saved in that file 'input'.
"""
import input
import matplotlib.pyplot as plt
import numpy as np

class DragForce:
    def __init__(self, v, t):
        self.t = t
        self.temperature = input.temperature
        self.atm_pressure = input.atm_pressure
        self.v = v # (m/s)
        self.radius = (input.a + input.b)/2

    def viscosity(self):
        return 2.791E-7 * self.temperature ** 0.7355

    def pressure(self):
        if input.pump == 'on':
            p = -10000 * self.t + self.atm_pressure
            if p > 200:
                return p
            elif p <= 200:
                return 200
        elif input.pump == 'off':
            return self.atm_pressure

    def air_density(self):
        R = 287.058
        return self.pressure()/(R*self.temperature)

    def Re(self):
        # calculate the Reynold number
        return self.air_density(self.t) * self.v * 2 * self.radius/self.viscosity()

    def mean_free_path(self):
         # Displacement of droplet moves between collisions
        l = np.sqrt(np.pi/8) * self.viscosity()/0.4987445 * 1/np.sqrt(self.air_density()*self.pressure())
        return l

    def cunningham_correction(self):
        Kn = self.mean_free_path()/self.radius
        alpha = 1.155
        beta = 0.471
        gamma = 0.596
        return 1 + Kn*(alpha + beta * np.exp(-gamma/Kn))

    def drag_force(self):
        F = 6 * np.pi * self.viscosity() * self.radius * self.v
        return F/self.cunningham_correction()

if __name__ == '__main__':
    # v = np.array([0.005,0.003])
    # DF = DragForce(300,v)
    # print(DF.atm_pressure)
    # print(DF.air_density())
    # # print(DF.cunningham_correction(0.011))
    # # print(DF.mean_free_path(0.011))
    # # print(DF.viscosity())
    # # print(DF.drag_force(0.011))

    # t = np.linspace(0,0.3,1000)
    # p = -10000 * t + 1000
    # for i in range(len(p)):
    #     if p[i] < 200:
    #         p[i] = 200
    # plt.figure('pressure vs time')
    # plt.plot(t, p)
    # plt.title('pressure vs time')
    # plt.xlabel('Time (s)')
    # plt.ylabel('pressure (Pa)')
    # plt.grid()
    # plt.show()
    pass