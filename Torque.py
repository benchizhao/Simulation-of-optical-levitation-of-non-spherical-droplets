# -*- coding: utf-8 -*-
"""
Time: 2020/Aug/26

@author: Benchi Zhao
Python version: 3.7
Functionality: This module is to calculate the torque caused by radiation pressure and air resistance.
Dependence: Need import 'Radiation_Force', 'drag_force',all parameters of the experiment is saved in that file 'input'.
"""

import numpy as np
import Radiation_Force as rf
import input
from drag_force import DragForce
def mass():
    return 4/3 * np.pi * input.a * input.b**2 * input.density

def momenrum_inertia():
    return 1/4 * mass() * (input.a**2 + input.b**2)

def rad_torque(rays, tilt, incident_angle):
    '''
    calculate the torque induced by radiation force
    :param rays: list
        all state of rays that interact with droplet
    :param tilt: list
        all tilt angle of interface
    :param incident_angle: list
        all incident angle between ray and interface
    :return: float
        total torque caused by the radiation force
    '''
    torque = 0
    for i in range(len(rays)):
        f1 = rf.rad_force_2_front(i=i,ray1=rays[i][5],ray2=rays[i][6],tilt=tilt[i][2],incident_angle=incident_angle[i][2])
        r1 = np.array(rays[i][5][:2]) - input.droplet_pos_e[:2]
        torque1 = np.cross(r1, f1)

        f2 = rf.rad_force_2_back(i=i,ray1=rays[i][7],ray2=rays[i][8],tilt=tilt[i][3],incident_angle=incident_angle[i][3])
        r2 = rays[i][7][:2] - input.droplet_pos_e[:2]
        torque2 = np.cross(r2, f2)
        torque = torque + torque1 + torque2
    return torque

def dra_torque(t):
    '''
    calculate the torque induced by air resistance
    :param t: float
        time which is used to solve ode
    :return:
        total torque caused by the drag force (air resistance)
    '''
    w = input.droplet_pos_e[5]
    r = np.sqrt(1/4*(input.a)**2 + 3/4*(input.b)**2)
    df = DragForce(v=w*r, t=t)
    f = abs(df.drag_force())
    if w > 0:
        return -r*f
    elif w <= 0:
        return r*f

def total_torque(rays, incident_angle, tilt, t):
    '''
    the total torque.
    :param rays: list
        all state of rays that interact with droplet
    :param tilt: list
        all tilt angle of interface
    :param incident_angle: list
        all incident angle between ray and interface
    :param t: float
        time which is used to solve ode
    :return: float
        total torque
    '''
    radiation_torque = rad_torque(rays, tilt, incident_angle)
    drag_torque = dra_torque(t)
    # print(radiation_torque, drag_torque)
    return radiation_torque + drag_torque

def angular_acc(rays, incident_angle, tilt, t):
    '''
    the angular acceleration.
    :param rays: list
        all state of rays that interact with droplet
    :param tilt: list
        all tilt angle of interface
    :param incident_angle: list
        all incident angle between ray and interface
    :param t: float
        time which is used to solve ode
    :return: float
        total torque
    '''
    return total_torque(rays, incident_angle, tilt, t)/momenrum_inertia()

if __name__ == '__main__':
    # print(total_torque())
    # print(np.shape(f_front))
    # print(np.shape(f_back))
    # print(mass())
    rays, incident_angle, tilt = rf.filter()
    t = total_torque(rays=rays, incident_angle=incident_angle, tilt=tilt, t=1)
    print(t)