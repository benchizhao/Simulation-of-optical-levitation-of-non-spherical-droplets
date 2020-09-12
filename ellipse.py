# -*- coding: utf-8 -*-
"""
Time: 2020/jul/20

@author: Benchi Zhao
Python version: 3.7
Functionality: This module is to find the interaction between straight line and ellipse, and to find the tangential line of the ellipse.
Dependence: Need import 'input', all parameters of the experiment is saved in that file 'input'.
"""
import numpy as np
import sympy
import matplotlib.pyplot as plt
import input
from scipy.optimize import fsolve
import timeit

def ellipse(a, b, center_point, theta):
    '''
    In Cartesian plane, any ellipse and write in the form of:
        Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
    :param a: float
        major axis of the ellipse
    :param b: float
        minor axis of the ellipse
    :param center_point: list,
        [x0,y0], the centre point of the ellipse, x0 and y0 are float
    :param theta: float
        the tilt angle, the angle from the position axis horizontal axis to the major axis
    :return: A, B, C, D, E, F floats
        These constant values describe the ellipse
    '''
    x0, y0 = center_point
    A = a**2 * np.sin(theta)**2 + b**2 * np.cos(theta)**2
    B = 2*(b**2 - a**2) * np.sin(theta) * np.cos(theta)
    C = a**2 * np.cos(theta)**2 + b**2 * np.sin(theta)**2
    D = -2*A*x0 - B*y0
    E = -B*x0 - 2*C*y0
    F = A*x0**2 + B*x0*y0 + C*y0**2 - (a*b)**2
    return A, B, C, D, E, F

def line(angle, point):
    '''
    line in Cartesian plane is y = kx + m
    :param angle: float
        angle between x-positive axis and the line
    :param point: list
        [x,y], one point on the line, x and y are float
    :return: k,b are the parameters to describe the line
    '''
    x, y = point
    k = np.tan(angle)
    m = y - k*x
    return k, m

def intersection(a, b, center_point, theta, angle, point):
    '''
    This function could find the intersection points between line and ellipse
    :param a: float
        major axis of the ellipse
    :param b: float
        minor axis of the ellipse
    :param center_point: list,
        [x0,y0], the centre point of the ellipse, x0 and y0 are float
    :param theta: float
        the tilt angle, the angle from the position axis horizontal axis to the major axis
    :param angle: float
        angle between x-positive axis and the line
    :param point: list
        [x,y], one point on the line, x and y are float
    :return: intersection points

    '''
    A, B, C, D, E, F = ellipse(a, b, center_point, theta)
    k, m = line(angle, point)
    # print(k,m)
    aa = A + B*k + C*k**2
    bb = B*m + 2*C*k*m + D + E*k
    cc = C*m**2 + E*m + F

    delta = bb**2 - (4 * aa*cc)

    if delta < 0:
        return []
    elif delta == 0:
        x = -bb/(2*aa)
        y = k * x + m
        return (x, y)
    else:
        x1 = (-bb - np.sqrt(delta))/(2*aa)
        y1 = k*x1 + m
        x2 = (-bb + np.sqrt(delta))/(2*aa)
        y2 = k * x2 + m
        return [[x1,y1], [x2,y2]]

def tangent_angle(intersection_pt):
    '''
    A straight forwards and slow method to find the ellipse tangential line which pass a given point, not used.
    :param intersection_pt: array
        a point array([x,y]) which is on the ellipse.
    :return: angle of the tangential line (angle between x-positive direction and the tangential line)
    '''
    x0, y0 = intersection_pt
    if y0 < 1E-10 and y0 > -1E-10:
        y0 = 0
    k = sympy.Symbol('k') #, real=True
    m = -k * x0 + y0
    A, B, C, D, E, F = ellipse(input.a, input.b, input.droplet_pos_e[:2], input.droplet_pos_e[2])
    AA = A + B * k + C * k ** 2
    BB = B * m + 2 * C * k * m + D + E * k
    CC = C * m ** 2 + E * m + F
    solved_value = sympy.solve([BB**2 - 4*AA*CC], [k])
    # print('solved_value', solved_value)
    if solved_value == []:
        return 0
    else:
        if str(solved_value[0][0])[0] == '-':
            value = -float(abs(solved_value[0][0]))
        else:
            value = float(abs(solved_value[0][0]))
        # print('value:', value)
        if value > 0:
            return np.pi/2 - abs(np.arctan(value))
        elif value < 0:
            return -np.pi/2 + abs(np.arctan(value))
        elif value ==0:
            return 0

def tangent_angle_by_derivation(intersection_pt):
    '''
    A complex and faster method to find the ellipse tangential line which pass a given point, which is used in the program.
    :param intersection_pt: array
        a point array([x,y]) which is on the ellipse.
    :return: angle of the tangential line (angle between x-positive direction and the tangential line)
    '''
    x0, y0 = intersection_pt
    if y0 < 1E-10 and y0 > -1E-10:
        y0 = 0
    A, B, C, D, E, F = ellipse(input.a, input.b, input.droplet_pos_e[:2], input.droplet_pos_e[2])
    if B*x0+2*C*y0+E == 0:
        return 0
    else:
        dfdx = (-D-2*A*x0-B*y0)/(B*x0+2*C*y0+E)
        slope = abs(np.pi/2 - abs(np.arctan(dfdx)))
        if dfdx >0:
            return -slope
        elif dfdx <0:
            return slope

if __name__ == '__main__':
    a = input.a
    b = input.b
    center_point = input.droplet_pos_e[:2]
    e_rot_angle = input.droplet_pos_e[2]

    l_angle = 0
    l_point = [900, 20E-6]
    inter = intersection(a, b, center_point, e_rot_angle, l_angle, l_point)
    tan = tangent_angle_by_derivation(inter[0])
    print(inter,'\n', tan)

    # print(tangent_angle_by_derivation(inter[0]))

