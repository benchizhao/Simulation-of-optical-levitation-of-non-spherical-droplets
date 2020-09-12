# -*- coding: utf-8 -*-
"""
Time: 2020/jul/10

@author: Benchi Zhao
Python version: 3.7
Functionality: This module is solve the trajectories of the droplet by ODE solvers
"""
import numpy as np
import input
import Radiation_Force as rf
import buoyancy as bf
import Torque
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import timeit
import pickle
from drag_force import DragForce

#prepare the forces
def mass():
    # mass of the droplet
    return 4/3 * np.pi * input.a * input.b**2 * input.density

def gravity():
    # prepare the gravity force for future use in diff_ellipse
    N = mass()*input.g
    return np.array([-N, float(0.)])

def drag_force(t):
    # prepare the drag force for future use in diff_ellipse
    v = input.droplet_pos_e[3:5]
    df = DragForce(v=v, t=t)
    drag_f = df.drag_force()
    return drag_f

def buoyancy(t):
    # prepare the buoyancy for future use in diff_ellipse
    v = input.droplet_pos_e[3:5]
    df = DragForce(v=v, t=t)
    air_density = df.air_density()
    buoyancy_f = bf.buoyancy(air_density)
    return buoyancy_f

def radiation_force(rays, incident_angle, tilt):
    # prepare the gravity force for future use in diff_ellipse
    return rf.radiation_force_2(rays, incident_angle, tilt)
    # return rf.radiation_force_1(rays)

def total_force(rays, incident_angle, tilt, t):
    # calculate the the total force that act on the droplet
    return gravity() + radiation_force(rays, incident_angle, tilt) - drag_force(t) + buoyancy(t)

def trans_acc(rays, incident_angle, tilt, t):
    # calculate the the acceleration of the droplet
    return total_force(rays, incident_angle, tilt, t)/mass()

def power_output(mode):
    if mode in ['c','const','constant']:
        return 2.5
    elif mode in ['noisy', 'n']:
        return 2.5 * (1 + (np.random.normal(loc=0.0, scale=1.0, size=None)) / 100000) # noise

#prepare the function for ode solver
def diff_ellipse(d_list, t):
    '''
    To solve the time ordinary derivative function.
    :param d_list: list
        state of the droplet [x,y,theta,vx,vy,angular_velocity]
    :param t: folat
        time
    :return: array
        d(d_list)/dt = [vx, vy, w, ax, ay, angular_acc]
    '''
    print('t=', t)
    input.droplet_pos_e = d_list

    input.power = power_output(input.power_mode)

    # if t > 0.1 and t < 0.101: # for square shape output
    #     input.power = 0.0000000001
    # elif t > 0.13 and t < 0.131:
    #     input.power = 2
    # else:
    #     input.power = 1

    # print('power:', input.power)
    rays, incident_angle, tilt = rf.filter()
    x, y, theta, vx, vy, w = d_list
    ax, ay = trans_acc(rays, incident_angle, tilt, t)
    angular_acc = Torque.angular_acc(rays, incident_angle, tilt, t)
    print('x=', x, 'y=', y, 'theta=', theta, 'vx=', vx, 'vy=', vy, 'w=', w)
    return np.array([vx, vy, w, ax, ay, angular_acc])


if __name__ == '__main__':
    start = timeit.default_timer()

    t = np.linspace(0, 0.15, 3000)

    result = odeint(diff_ellipse, input.droplet_pos_e, t)
    f = open('pump air out-2.txt', 'wb')
    pickle.dump(result, f)
    f.close()



    '''
    Plot 4 figures: x-position, x-velocity, y-position, y-velocity respect to time
    '''
    f = open('pump air out-2.txt', 'rb')
    d = pickle.load(f)
    f.close()

    f = open('circle all on new code.txt', 'rb')
    e = pickle.load(f)
    f.close()

    # print(d)

    # t = np.linspace(0, 0.05, 3000)

    plt.figure('')
    plt.plot(t, d[:, 0],'-', label='pump air out')
    plt.plot(t, e[:, 0],'-.',label='constant air pressure' )

    # plt.title('x,y,theta coupling')
    plt.xlabel('Time (s)')
    plt.ylabel('x')
    # plt.ylim(0,900E-6)
    plt.grid()
    plt.legend()
    plt.show()

    #
    # plt.figure('x-position of droplet power_off[0.15,0.16]')
    # plt.plot(t, d[:, 0])
    # plt.title('x-position of droplet')
    # plt.xlabel('Time (s)')
    # plt.ylabel('x-axis')
    # # plt.ylim(0,900E-6)
    # plt.grid()
    # plt.show()
    #
    # plt.figure('y-position of droplet power_off[0.15,0.16]')
    # plt.plot(t, d[:, 1])
    # plt.title('y-position of droplet')
    # plt.xlabel('Time (s)')
    # plt.ylabel('y-axis')
    # plt.grid()
    # plt.show()
    #
    # plt.figure('x_velocity of droplet power_off[0.15,0.16]')
    # plt.plot(t, d[:, 2])
    # plt.title('x_velocity of droplet')
    # plt.xlabel('Time (s)')
    # plt.ylabel('x-velocity')
    # plt.grid()
    # plt.show()
    #
    # plt.figure('y_velocity of droplet power_off[0.15,0.16]')
    # plt.plot(t, d[:, 3])
    # plt.title('y_velocity of droplet')
    # plt.xlabel('Time (s)')
    # plt.ylabel('y-velocity')
    # plt.grid()
    # plt.show()

    # plt.figure('x position respect to x velocity')
    # plt.plot(d[:,0], d[:, 2])
    # plt.title('x position respect to x velocity')
    # plt.xlabel('x position')
    # plt.ylabel('x velocity')
    # plt.grid()
    # plt.show()

    # plt.figure('y position respect to y velocity')
    # plt.plot(d[:, 1], d[:, 3])
    # plt.title('y position respect to y velocity')
    # plt.xlabel('y position')
    # plt.ylabel('y velocity')
    # plt.grid()
    # plt.show()

    stop = timeit.default_timer()
    print('Time: ', stop - start)

