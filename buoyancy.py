# -*- coding: utf-8 -*-
"""
Time: 2020/Aug/20

@author: Benchi Zhao
Python version: 3.7
Functionality: This module calculate the buoyancy force in the system, which is a constant force encounters the gravity.
Dependence: Need import 'input', all parameters of the experiment is saved in that file 'input'.
"""

import numpy as np
import input

def volume():
    return 4/3*np.pi * input.a * input.b**2
def buoyancy(air_density):
    '''
    calculate the buoyancy force.
    :param air_density: float
        the air density of the surrounding
    :return:
    '''
    f = air_density * input.g * volume()
    return np.array([f, 0])

print(volume())