# -*- coding: utf-8 -*-
"""
Time: 2020/Jul/06

@author: Benchi Zhao
Python version: 3.7
Functionality: This module trace one individual ray and plot the diagram.
Dependence: Need import 'input', all parameters of the experiment is saved in that file 'input'.
"""
import input
import ellipse as ep
import pandas as pd
import numpy as np
import math
import copy

import matplotlib.pyplot as plt
import timeit
import sys
from numba import jit
import timeit
from multiprocessing import Pool

def degree2rad(degree):
    return degree/180 * math.pi
def rad2degree(rad):
    return rad/math.pi * 180
def slope2rad(slope):
    return np.arctan(slope)
def rad2slope(rad):
    return np.tan(rad)
def degree2slope(degree):
    return np.tan(np.deg2rad(degree))
def slope2degree(slope):
    return np.rad2deg(np.arctan(slope))

def gaussian(y, sigma):
    '''
    :param y: float
        the y position describes the distance between beam axis and the the ray
    :return: float
        the gaussian value at y
    '''
    return 1 / (np.sqrt(2 * math.pi)*sigma) * math.exp(-y ** 2 / (2*sigma**2))
def snell(theta1):
    '''
           snells(theta1)
               Calculates theta2 from snells law in radians.

           Parameters
           ----------
           theta1: float
               Angle in to refracting surface (in rad).

           Returns
           -------
           theta2: float
               Angle out of refracting surface.

           inputut Script Parameters
           ---------------------------
           medium_n: float
               Refractive index of the surrounding medium.
           target_n: float
               Refractive index of target.

           '''
    theta2 = np.arcsin((input.medium_n / input.target_n) * np.sin(theta1))
    return theta2

def R_value(theta1): #in rad
    '''
    R(theta1)
        Calculates the Fresnel reflectivity as a function of theta1.

    Parameters
    -------------
    theta1: float
        Angle in to refracting surface.

    Returns
    ---------
    R: float
        Fresnel power reflectance as a function of theta1 and theta2.

    Input Script Parameters
    ---------------------------
    polarisation: string
        Polarisation of laser beam used.
    medium_n: float
        Refractive index of the surrounding medium.
    target_n: float
        Refractive index of target.

    '''
    theta2 = snell(theta1)
    if input.polarisation == 'p':
        result = (abs((input.medium_n * np.cos(theta2) - input.target_n * np.cos(theta1)) / (
                    input.medium_n * np.cos(theta2) + input.target_n * np.cos(theta1)))) ** 2
        return result
    elif input.polarisation == 's':
        result = ((input.medium_n*np.cos(theta1)-(input.target_n*np.cos(theta2)))/(input.medium_n*np.cos(theta1)+(input.target_n*np.cos(theta2))))**2
        return result
    elif input.polarisation == 'circle':
        result = 0.5*(((np.sin(theta1-theta2)**2)/(np.sin(theta1+theta2)**2))+((np.tan(theta1-theta2)**2)/(np.tan(theta1+theta2)**2)))
        return result
    else:
        sys.exit("Polarisation is not defined it must be set to 'p','circle' or 's', check input file")
def T_value(theta1):
    '''
    T(theta1)
        Calculates the Fresnel reflectivity as a function of theta1.

    Parameters
    -------------
    theta1: float
        Angle in to refracting surface.

    Returns
    ---------
    R: float
        Fresnel refractivity as a function of theta1.

    Input Script Parameters
    ---------------------------
    polarisation: string
        Polarisation of laser beam used.
    medium_n: float
        Refractive index of the surrounding medium.
    target_n: float
        Refractive index of target.

    '''
    if input.target_n >= input.medium_n:
        reflection = R_value(theta1)
    else:
        thetac = np.arcsin(input.target_n / input.medium_n)
        if theta1 >= thetac:
            reflection = 1
        else:
            reflection = R_value(theta1)
    T = 1 - reflection
    return T

def circle_line_segment_intersection(circle_center, circle_radius, pt1, pt2, full_line=True, tangent_tol=1e-9):
    """ Find the points at which a circle intersects a line-segment.  This can happen at 0, 1, or 2 points.

    :param circle_center: The (x, y) location of the circle center
    :param circle_radius: The radius of the circle
    :param pt1: The (x, y) location of the first point of the segment
    :param pt2: The (x, y) location of the second point of the segment
    :param full_line: True to find intersections along full line - not just in the segment.  False will just return intersections within the segment.
    :param tangent_tol: Numerical tolerance at which we decide the intersections are close enough to consider it a tangent
    :return Sequence[Tuple[float, float]]: A list of length 0, 1, or 2, where each element is a point at which the circle intercepts a line segment.

    Note: We follow: http://mathworld.wolfram.com/Circle-LineIntersection.html
    """

    (p1x, p1y), (p2x, p2y), (cx, cy) = pt1, pt2, circle_center
    (x1, y1), (x2, y2) = (p1x - cx, p1y - cy), (p2x - cx, p2y - cy)
    dx, dy = (x2 - x1), (y2 - y1)
    dr = (dx ** 2 + dy ** 2)**.5
    big_d = x1 * y2 - x2 * y1
    discriminant = circle_radius ** 2 * dr ** 2 - big_d ** 2

    if discriminant < 0:  # No intersection between circle and line
        return []
    else:  # There may be 0, 1, or 2 intersections with the segment
        intersections = [
            (cx + (big_d * dy + sign * (-1 if dy < 0 else 1) * dx * discriminant**.5) / dr ** 2,
             cy + (-big_d * dx + sign * abs(dy) * discriminant**.5) / dr ** 2)
            for sign in ((1, -1) if dy < 0 else (-1, 1))]  # This makes sure the order along the segment is correct
        if not full_line:  # If only considering the segment, filter out intersections that do not fall within the segment
            fraction_along_segment = [(xi - p1x) / dx if abs(dx) > abs(dy) else (yi - p1y) / dy for xi, yi in intersections]
            intersections = [pt for pt, frac in zip(intersections, fraction_along_segment) if 0 <= frac <= 1]
        if len(intersections) == 2 and abs(discriminant) <= tangent_tol:  # If line is tangent to circle, return just one point (as both intersections have same location)
            return [intersections[0]]
        else:
            return intersections

class Ray():
    '''
    Simulate the path of rays
    '''
    def __init__(self):
        self.droplet_position = input.droplet_pos_e
        self.ray_state = []
        self.incident_angle = []
        self.tilt = []

    def trigger(self, y, intensity):
        '''
            Originate the beam at the beginning
            4 parameters are used to describe the ray [x,y,theta,intensity];

        x: float
            x position of the initial beam.
        y: float
            y position of the initial beam.
        theta: float
            the tilt angle of the beam. (in rad)
        intensity: float
            The intensity of the beam.
        '''
        if input.ray_type in ['Gaussian', 'gaussian', 'G', 'g']:
            ray = [0, y, 0, intensity]
            self.ray_state.append(ray)
        elif input.ray_type in ['Flat', 'flat', 'f', 'F']:
            ray = [0, y, 0, 1]
            self.ray_state.append(ray)
        else:
            raise NameError

    def free_propagate(self, distance):
        '''
        Describe the ray state after the ray propagate for a distance
        :param distance: float
            The distance that the ray propagate
        :return list
            [x,y,theta,intensity]
        '''
        state = self.ray_state[-1]
        # if state[2] < 0:
        #     self.ray_state[-1][2] += 2 * np.pi

        if state[2] == np.pi / 2:
            ray = [state[0], state[1] + distance, state[2], state[3]]
        elif state[2] == 3 * np.pi / 2:
            ray = [state[0], state[1] - distance, state[2], state[3]]
        elif state[2] > np.pi / 2 and state[2] < 3 * np.pi / 2:
            slope = rad2slope(state[2])
            x = state[0] - distance
            y = -slope * distance + state[1]
            ray = [x, y, state[2], state[3]]
        else:
            slope = rad2slope(state[2])
            x = state[0] + distance
            y = slope * distance + state[1]
            ray = [x, y, state[2], state[3]]
        return ray

    def full_reflection(self, tilt):
        '''
        append the ray state after full reflection
        :param tilt: float
            The angle tilted by the interface, if the tile is 0, then the interface is verticle;
                                               if the tile is positive, the interface rotate anti-clockwise;
                                               if the tilt is negative, the interface rotate clockwise. (in rad)
        :return:
            ray [x,y,theta,intensity]
        '''
        k = np.tan(np.pi/2 + tilt) # slope of the mirror
        x0, y0 = self.ray_state[-1][:2] # point on the inter surface
        a, b = self.ray_state[-2][:2]
        A = -k
        B = 1
        C = k*x0 - y0
        symmetry_pt = [-(2*A*B*b+(A**2-B**2)*a+2*A*C)/(B**2+A**2), -((B**2-A**2)*b+2*A*B*a + 2*B*C)/(B**2+A**2)]
        k_new = (symmetry_pt[1]-y0)/(symmetry_pt[0]-x0)
        # print('check:', k*(symmetry_pt[1]-b)/(symmetry_pt[0]-a))
        if symmetry_pt[1] > y0:
            angle = -abs(np.arctan(k_new))
        else:
            angle = np.pi - abs(np.arctan(k_new))
        ray = [self.ray_state[-1][0], self.ray_state[-1][1], angle, self.ray_state[-1][3]]
        return ray

    def flat_interface(self, tilt, n_left, n_right):
        '''
        Describe how the ray state change after the ray pass through the flat interface.
        Append the refractive beam into forward, and append the reflected beam into backward.

        parameter
        -------------
        tilt: float
            The angle tilted by the interface, if the tile is 0, then the interface is verticle;
                                               if the tile is positive, the interface rotate anti-clockwise;
                                               if the tilt is negative, the interface rotate clockwise. (in rad)
        n_left: float
            The refractive index of the incident beam.
        n_right: float
            The refractive index of the refractive beam.
        -----------------------------------------------
        :return list
            The ray state after going through the flat interface [x, y, theta, intensity]
        '''
        state = self.ray_state[-1]
        # snell's law method
        incident_rad = state[2] - tilt
        self.tilt.append(tilt)
        self.incident_angle.append(incident_rad)
        if abs(n_left * np.sin(incident_rad)/n_right) >= 1: # full reflection
            ray = self.full_reflection(tilt)
            if ray[2] >=np.pi:
                ray[2] = ray[2] - np.pi
            # ray = [state[0], state[1], np.pi - state[2], state[3]]
        else:
            refracted_deg = np.arcsin(n_left * np.sin(incident_rad)/n_right) + tilt
            # reflected_angle = np.pi - state[2] + 2 * tilt
            if refracted_deg>=np.pi:
                refracted_deg = refracted_deg- np.pi
            ray = [state[0], state[1], refracted_deg, state[3] * T_value(incident_rad)]
        return ray

    def lens(self, f, position, thick):
        '''
        The ray goes through the two inferfaces of the lens (into the droplet and out of the droplet).
        Describe the ray state after going through the lens interfaces, and two ray states are appended into self.ray_state list
        :param f: float
            focal length of the droplet
        :param position: float
            central point position of the lens
        :param thick:
            thickness(width) of the lens
        -----------------------------------------------
        :return:
            Append the two ray states after lens into the self.ray_state list
        '''
        r = 2*f*(input.target_n-input.medium_n)

        self.ray_state.append(self.free_propagate(position))
        center1 = [position-(r-thick/2), 0]
        center2 = [position + (r - thick / 2), 0]
        pt1 = copy.deepcopy(self.ray_state[-1][0:2])
        pt2 = [pt1[0]+1, pt1[1]+np.tan(self.ray_state[-1][2])]
        left_intersection_p = circle_line_segment_intersection(center2, r, pt1, pt2, full_line=True, tangent_tol=1e-9)[0]
        self.ray_state[-1][:2] = left_intersection_p # hit the left surface

        tilt = np.arcsin((self.ray_state[-1][1] - center1[1])/ r)
        self.ray_state.append(self.flat_interface(tilt=-tilt, n_left=input.medium_n, n_right=input.target_n)) # begin to propagate in the lens

        pt1 = copy.deepcopy(self.ray_state[-1][0:2])
        pt2 = [pt1[0]+1, pt1[1]+np.tan(self.ray_state[-1][2])]
        right_intersection_p = circle_line_segment_intersection(center1, r, pt1, pt2, full_line=True, tangent_tol=1e-9)[1]
        self.ray_state.append([right_intersection_p[0], right_intersection_p[1], self.ray_state[-1][2], self.ray_state[-1][3]])

        tilt = np.arcsin((self.ray_state[-1][1] - center1[1]) / r)
        self.ray_state.append(self.flat_interface(tilt=tilt, n_left=input.target_n, n_right=input.medium_n))

    def circle_droplet(self, radius):
        '''
        The ray goes through the two inferfaces of the circle droplet (into the droplet and out of the droplet).
        Describe the ray state after going through the circle droplet interfaces, and two ray states are appended into self.ray_state list
        :param radius: float
            radius of the droplet
        -----------------------------------------------
        :return:
            Append two ray states into the self.ray_state list
        '''
        position = input.droplet_pos_e[:2]
        self.ray_state.append(self.free_propagate(position[0]-input.lens_pos))
        pt1 = [self.ray_state[-1][0], self.ray_state[-1][1]]
        pt2 = [pt1[0]+1, pt1[1]+np.tan(self.ray_state[-1][2])]
        intersection = circle_line_segment_intersection(position, radius, pt1, pt2, full_line=True, tangent_tol=1e-13)

        if intersection == []:
            self.free_propagate(radius)
        else:
            intersection = intersection[0]
            self.ray_state[-1][:2] = intersection
            tilt = np.arcsin((self.ray_state[-1][1] - input.droplet_pos_e[1]) / radius)
            self.ray_state.append(self.flat_interface(tilt=-tilt, n_left=input.medium_n, n_right=input.target_n))

            pt1 = [self.ray_state[-1][0], self.ray_state[-1][1]]
            pt2 = [pt1[0] + 1, pt1[1] + np.tan(self.ray_state[-1][2])]
            intersection = circle_line_segment_intersection(position, radius, pt1, pt2, full_line=True, tangent_tol=1e-13)[1]
            self.ray_state.append([intersection[0], intersection[1], self.ray_state[-1][2], self.ray_state[-1][3]])
            tilt = np.arcsin((self.ray_state[-1][1] - input.droplet_pos_e[1]) / radius)
            self.ray_state.append(self.flat_interface(tilt=tilt, n_left=input.target_n, n_right=input.medium_n))
            # self.ray_state.append(self.flat_interface(tilt=tilt, n_left=input.medium_n, n_right=input.target_n))

        self.ray_state.append(self.free_propagate(100E-6))

    def ellipse_droplet(self, a, b, center_point, theta):
        '''
        The ray goes through the two inferfaces of the ellipse droplet (into the droplet and out of the droplet).
        Describe the ray state after going through the ellipse droplet interfaces, and two ray states are appended into self.ray_state list
        :param a: float
            major axis of the ellipse
        :param b: float
            minor axis of the ellipse
        :param center_point: list,
            [x0,y0], the centre point of the ellipse, x0 and y0 are float
        :param theta: float
            the tilt angle, the angle from the position axis horizontal axis to the major axis
        -----------------------------------------------
        :return:
            Append two ray states into the self.ray_state list
        '''
        # light propagate from lens to droplet
        position = input.droplet_pos_e[:2]
        self.ray_state.append(self.free_propagate(position[0] - input.lens_pos))

        # incident into the droplet
        x, y, angle = self.ray_state[-1][:3]
        inter_point1 = ep.intersection(a, b, center_point, theta, angle, [x, y])



        if len(inter_point1) != 2: # not interact with droplet
            self.ray_state.append(self.free_propagate(a))
        else: # interact with droplet
            inter_point1 = inter_point1[0]
            self.ray_state[-1][:2] = inter_point1

            tilt = ep.tangent_angle_by_derivation(inter_point1)
            self.ray_state.append(self.flat_interface(tilt=tilt, n_left=input.medium_n, n_right=input.target_n))
            inter_point2 = ep.intersection(a, b, center_point, theta, self.ray_state[-1][2], self.ray_state[-1][:2])[1]
            self.ray_state.append(self.free_propagate(inter_point2[0] - inter_point1[0]))
            tilt = ep.tangent_angle_by_derivation(inter_point2)
            self.ray_state.append(self.flat_interface(tilt=tilt, n_left=input.target_n, n_right=input.medium_n))
            # self.ray_state.append(self.flat_interface(tilt=-tilt, n_left=input.medium_n, n_right=input.target_n))

        self.ray_state.append(self.free_propagate(100E-6))

    def make_table(self):
        '''
        Show the all_forwards and all_backwards in table form,  making it clean.
        '''
        # for i in range(4):
        #     data = self.ray_state[i]
        #     print(np.shape(data))
        df_1 = pd.DataFrame(self.ray_state, columns=['x', 'y', 'theta', 'intensity'])
        print(df_1)

    def plot_single_ray_c(self):
        '''
        Plot the path of single ray that interacting with lens and circular droplet
        '''
        plt.figure(1)
        for i in range(len(self.ray_state)-1):
            x = np.linspace(self.ray_state[i][0], self.ray_state[i+1][0], 10)
            y = np.linspace(self.ray_state[i][1], self.ray_state[i+1][1], 10)
            plt.plot(x, y, 'k')
        # plot lens
        theta1 = np.linspace(3 * math.pi / 4, 5 * math.pi / 4)
        theta2 = np.linspace(- math.pi / 4, math.pi / 4)
        x_1 = input.lens_f * np.cos(theta1) + input.lens_pos + input.lens_f - input.len_thickness / 2
        x_2 = input.lens_f * np.cos(theta2) + input.lens_pos - input.lens_f + input.len_thickness / 2
        y_1 = input.lens_f * np.sin(theta1)
        y_2 = input.lens_f * np.sin(theta2)
        plt.plot(x_1, y_1, 'g')
        plt.plot(x_2, y_2, 'g')
        #plot circle droplet
        theta = np.linspace(0,2*np.pi)
        x = input.radius * np.cos(theta) + input.droplet_pos_e[0]
        y = input.radius * np.sin(theta) + input.droplet_pos_e[1]
        plt.plot(x,y)
        plt.ylim(-50E-6, 50E-6)
        plt.show()

    def plot_single_ray_e(self):
        '''
        Plot the path of single ray that interating with lens and ellipse droplet
        '''
        plt.figure(1)
        for i in range(len(self.ray_state) - 1):
            x = np.linspace(self.ray_state[i][0], self.ray_state[i + 1][0], 10)
            y = np.linspace(self.ray_state[i][1], self.ray_state[i + 1][1], 10)
            plt.plot(x, y, 'k')
        # plot lens
        theta1 = np.linspace(3 * math.pi / 4, 5 * math.pi / 4)
        theta2 = np.linspace(- math.pi / 4, math.pi / 4)
        x_1 = input.lens_f * np.cos(theta1) + input.lens_pos + input.lens_f - input.len_thickness / 2
        x_2 = input.lens_f * np.cos(theta2) + input.lens_pos - input.lens_f + input.len_thickness / 2
        y_1 = input.lens_f * np.sin(theta1)
        y_2 = input.lens_f * np.sin(theta2)
        plt.plot(x_1, y_1, 'grey')
        plt.plot(x_2, y_2, 'grey')
        # plot ellipse droplet
        theta = np.linspace(0, 2 * np.pi, 100)
        a=input.a
        b=input.b
        e_rot_angle = input.droplet_pos_e[2]
        center_point = input.droplet_pos_e[:2]
        x = a * np.cos(theta)
        y = b * np.sin(theta)
        #
        x_e = np.cos(e_rot_angle) * x - np.sin(e_rot_angle) * y + center_point[0]
        y_e = np.sin(e_rot_angle) * x + np.cos(e_rot_angle) * y + center_point[1]
        plt.plot(x_e, y_e)
        plt.ylim(-50E-6, 50E-6)
        plt.show()

if __name__ == '__main__':
    '''Trace one individual ray'''

    start = timeit.default_timer()
    r = Ray()
    r.trigger(y=-1e-06, intensity=1) # the initial position of the ray
    r.lens(f=input.lens_f, position=input.lens_pos, thick=input.len_thickness) #the ray passes the convex lens

    if input.droplet_shape in ['c', 'circle']:  #the ray passes the circle
        r.circle_droplet(input.radius)
        print(r.ray_state)
        r.plot_single_ray_c()
    elif input.droplet_shape in ['e', 'ellipse']:
        r.ellipse_droplet(a=input.a, b=input.b, center_point=input.droplet_pos_e[:2], theta=input.droplet_pos_e[2]) #the ray passes the droplet
        print(r.ray_state)
        r.plot_single_ray_e()
    # r.make_table()


    # print(r.incident_angle)
    # print(r.tilt)
    end = timeit.default_timer()
    print('time:', end-start)
   