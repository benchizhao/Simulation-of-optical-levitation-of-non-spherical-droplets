# -*- coding: utf-8 -*-
"""
Time: 2020/Jul/06

@author: Benchi Zhao
Python version: 3.7
Functionality: This module determines the parameters used in the simulation.
"""

import numpy as np
'''
Ray Tracing parameters
----------------------
'''
# About the beam
ray_type = 'g'  # energy distribution of rays. 'f' represents flat, the intensity of all rays are the same;
                                            #  'g' represent gaussian, the intensity distriburion of rays are gaussian.
polarisation = 'p'  # Type of polarisation: 's', 'p', 'circle'
power = 1000E-3 # (W)Total power of the beam.
power_mode = 'c' # 'c' or 'constant' for constant power output; 'n' or 'noisy' for noisy output
sigma = 10E-6  # sigma is the parameter to describe the gaussian.
no_of_rays = 400   # Number of rays we are going to trace
width = 40E-6  # Width of the rays
# Refractive index
medium_n = 1    # Refractive index of the surrounding
target_n = 1.5  # Refractive index of the droplet
y_displacement = 0


'''
Parameters for droplet
'''
droplet_shape = 'c' # shape of the droplet, options are circle(c) and ellipse(e)
density = 960 # kg/m^3


# circle
radius = 20E-6   # (um) radius of the droplet
droplet_pos_c = np.array([950, 0, 0, 0])*1E-6      # (um)Central position of the circle droplet [x, y, vx, vy]


# ellipse
droplet_pos_e = np.array([1000*1E-6, 20*1E-6, 0, 0, 0, 0])      # (um)Central position of the ellipse droplet [x, y, tilt_angle, vx, vy, angular_v]
a = 20E-6 # (um) major axis for ellipse
b = 20E-6 # (um) minor axis for ellipse



'''
Parameters for lens
'''
lens_pos = 50E-6               # (um)Central position of the lens
lens_f = 200E-6                # (um)Focal length of the lens
len_thickness = 20E-6          # (um)Thickness of the lens

'''
Viscosity drag
'''
rho = 1.293     # Air density   kg/m^3
n = 1.8E-5    # Air viscosity kg/(m.s)

'''
Other useful parameters
'''
c = 2.98E+8 # Speed of light
g = 9.8
temperature = 300 # (K)temeprature of lab
atm_pressure = 101325 # (Pa)pressure of the environment

pump = 'on' # if the pump is on, the air pressure decay respects to time; if is off, the air pressure is constant
