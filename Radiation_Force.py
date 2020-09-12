# -*- coding: utf-8 -*-
"""
Time: 2020/Jul/26

@author: Benchi Zhao
Python version: 3.7
Functionality: This module calculate the force caused by radiation pressure.
Dependence: Need import 'NewRT', all parameters of the experiment is saved in that file 'input'.
"""

import numpy as np
import input
import RaysTracing as RsT
import RayTracing as RT
import timeit
import matplotlib.pyplot as plt
import Torque



def filter():
    '''
    Delete the rays that do not interact with droplet, and only the rays that interact with droplet remains.
    :return: list
        rays that interact with droplet
    '''
    rays, incident_angle, tilt = RsT.bundle() # achieves the state of rays
    for i in range(input.no_of_rays-1, -1, -1):
        if len(rays[i]) < 10:
            rays.remove(rays[i])
    for i in range(input.no_of_rays-1, -1, -1):
        if len(incident_angle[i]) < 4:
            incident_angle.remove(incident_angle[i])
    for i in range(input.no_of_rays-1, -1, -1):
        if len(tilt[i]) < 4:
            tilt.remove(tilt[i])
    return rays, incident_angle, tilt
'''
First method to calculate the forces on droplet (Ashkin):
This method is introduced by Ashkin's paper 'Forces of a Single-Beam Gradient Laser Trap on a Dielectric Sphere in the Ray Optics'
This method works only for 
'''
def F_s(data):
    '''
    Calculate the scattering force.
    :param data: list
        the data that the rays state interacting with droplet.
    :return total_force: list
        The total_foce contains two element, the first one is the force along x-aixs, and the second one is the force along y-axis.
    '''

    total_force = np.array([0, 0])
    for i in range(len(data)):
        x = data[i][5][0]
        y = data[i][5][1]
        vec_1 = [x-input.droplet_pos_e[0], y-input.droplet_pos_e[1]]
        vec_2 = [-1, -np.tan(data[i][5][2])]
        unit_vector_1 = vec_1 / np.linalg.norm(vec_1)
        unit_vector_2 = vec_2 / np.linalg.norm(vec_2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        theta1 = np.arccos(dot_product)
        theta2 = RT.snell(theta1)  # Refracted angle

        P = input.power * data[i][5][3]
        F = input.medium_n * P / input.c*(1+RT.R_value(theta1)*np.cos(2*theta1)- (RT.T_value(theta1) ** 2 * (np.cos(2 * theta1 - 2 * theta2) + RT.R_value(theta1) * np.cos(2 * theta1))) / (1 + RT.R_value(theta1) ** 2 + 2 * RT.R_value(theta1) * np.cos(2 * theta2)))
        Fx = abs(F) * np.cos(data[i][5][2])
        Fy = abs(F) * np.sin(data[i][5][2])
        total_force = total_force + np.array([Fx, Fy])
    return total_force

def F_g(data):
    '''
    Calculate the gradient force.
    :param data: list
        the data that the rays state interacting with droplet.
    :return total_forcr: list
        The total_foce contains two element, the first one is the force along x-aixs, and the second one is the force along y-axis.
    '''
    total_force = np.array([0, 0])
    for i in range(len(data)):
        x = data[i][5][0]
        y = data[i][5][1]
        vec_1 = [x-input.droplet_pos_e[0], y-input.droplet_pos_e[1]]
        vec_2 = [-1, -np.tan(data[i][5][2])]
        unit_vector_1 = vec_1 / np.linalg.norm(vec_1)
        unit_vector_2 = vec_2 / np.linalg.norm(vec_2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        theta1 = np.arccos(dot_product)
        theta2 = RT.snell(theta1) # Refracted angle

        P = input.power * data[i][5][3]
        F = input.medium_n*P/input.c *(RT.R_value(theta1)*np.sin(2*theta1)-(RT.T_value(theta1)**2*(np.sin(2*theta1-2*theta2)+RT.R_value(theta1)*np.sin(2*theta1)))/(1+RT.R_value(theta1)**2+2*RT.R_value(theta1)*np.cos(2*theta2)))
        if data[i][5][2] < np.pi and data[i][5][2]>0:
            Fx = abs(F) * np.cos(data[i][5][2]-np.pi/2)
            Fy = abs(F) * np.sin(data[i][5][2]-np.pi/2)
        else:
            Fx = abs(F) * np.cos(data[i][5][2] + np.pi / 2)
            Fy = abs(F) * np.sin(data[i][5][2] + np.pi / 2)
        total_force = total_force + np.array([Fx, Fy])
    return total_force

def radiation_force_1(data):
    '''
    Calculate the total force.
    :param data: list
        the data that the rays state interacting with droplet.
    :return total_forcr: list
        The total_force calculated by Ashkin's method.
    '''
    Force = F_g(data) + F_s(data)
    return Force

'''
Second method to calculate the forces on droplet:
The  equation used is F = P/c, where F is the force, P is the power of rays, c is the speed of light
'''

def rad_force_2_front(i, ray1, ray2, tilt, incident_angle):
    '''
    There are two interaction points between rays and droplets (into & out of), the function calculate the force
    induced by the i_th ray interacted on the first interaction point.
    :param i: int
        the index of rays, or which ray that is calculated.
    :param ray1: list
        the ray state of the i_th incident ray.
    :param ray2: list
        the ray state of the i_th refracted ray.
    :param tilt: float
        the angle of the interface tilt (rad)
    :param incident_angle: float
        the angle of the ray tilt (rad)

    :return: list
        The force induced by the i_th ray at the first interaction point.
    '''
    c = input.c
    power1 = ray1[3] * input.power
    power2 = ray2[3] * input.power
    power3 = power1 - power2

    F1 =np.array([power1 / c * np.cos(incident_angle), power1 / c * np.sin(incident_angle)])
    refract_angle = np.arcsin(input.medium_n/input.target_n*np.sin(incident_angle))
    F2 = np.array([power2 / c * np.cos(refract_angle + np.pi), power2 / c * np.sin(refract_angle+ np.pi)])
    F3 = np.array([F1[0] * power3/power1, -F1[1] * power3/power1])

    Fs, Fg = F1 + F2 + F3
    # rotation_matrix = np.mat([[np.cos(tilt), -np.sin(tilt)], [np.sin(tilt), np.cos(tilt)]])
    rotation_matrix = np.mat([[np.cos(tilt), np.sin(tilt)], [-np.sin(tilt), np.cos(tilt)]])
    force = rotation_matrix*np.mat([[Fs], [Fg]])
    force = force.flatten().A
    Fx = force[0][0]
    Fy = force[0][1]
    return np.array([Fx,Fy])

def rad_force_2_back(i, ray1, ray2, tilt, incident_angle):
    '''
    There are two interaction points between rays and droplets (into & out of), the function calculate the force
    induced by the i_th ray interacted on the second interaction point.
    :param i: int
        the index of rays, or which ray that is calculated.
    :param ray1: list
        the ray state of the i_th incident ray.
    :param ray2: list
        the ray state of the i_th refracted ray.
    :param tilt: float
        the angle of the interface tilt (rad)
    :param incident_angle: float
        the angle of the ray tilt (rad)

    :return: list
        The force induced by the i_th ray at the second interaction point.
    '''

    c = input.c
    power1 = ray1[3] * input.power
    power2 = ray2[3] * input.power
    power3 = power1 - power2

    F1 =np.array([power1 / c * np.cos(incident_angle), power1 / c * np.sin(incident_angle)])
    F2 = np.array([0, 0])
    if power1 == power2: #full reflection
        F3 = np.array([F1[0] * power2 / power1, -F1[1] * power2 / power1])
    else:
        refract_angle = np.arcsin(input.target_n/input.medium_n*np.sin(incident_angle))
        F2 = np.array([power2 / c * np.cos(refract_angle + np.pi), power2 / c * np.sin(refract_angle+ np.pi)])
        F3 = np.array([F1[0] * power3 / power1, -F1[1] * power3 / power1])
    Fs, Fg = F1 + F2 + F3
    # rotation_matrix = np.mat([[np.cos(tilt), -np.sin(tilt)], [np.sin(tilt), np.cos(tilt)]])
    rotation_matrix = np.mat([[np.cos(tilt), np.sin(tilt)], [-np.sin(tilt), np.cos(tilt)]])
    force = rotation_matrix * np.mat([[Fs], [Fg]])
    force = force.flatten().A
    Fx = force[0][0]
    Fy = force[0][1]
    return np.array([Fx, Fy])

def radiation_force_2(rays, incident_angle, tilt):
    '''
    The total force induced by the bundle of rays calculated by the second method
    :param rays: list
        all of the incident ray states.
    :param incident_angle: float
        all of the angles of the ray tilt (rad)
    :param tilt: float
        all of the angle of the interface tilt (rad)

    :return: list
        The total force induced by a bundle of rays by the second method.
    '''
    force = np.array([0,0])
    for i in range(len(rays)):
        f_front = rad_force_2_front(i=i,ray1=rays[i][5],ray2=rays[i][6],tilt=tilt[i][2],incident_angle=incident_angle[i][2])
        f_back = rad_force_2_back(i=i,ray1=rays[i][7],ray2=rays[i][8],tilt=tilt[i][3],incident_angle=incident_angle[i][3])
        force = force + f_front + f_back
        # print(i, f_front, f_back)
    # print('tot',force)
    return force


if __name__ == '__main__':
    def inter1(x_range, y_range): # plot the force (calculated by Ashkin's method) respect to position
        x = np.linspace(x_range[0], x_range[1], 300)
        y = np.linspace(y_range[0], y_range[1], 300)
        plt.figure()
        for i in range(len(x)):
            input.droplet_pos_e = [x[i], y[i], 0, 0, 0, 0]
            data, incident_angle, tilt = filter()
            print(i)
            if x_range[0] == x_range[1]:
                plt.plot(y[i], radiation_force_1(data)[1], 'k.')
                plt.title('y position respect to force')
                plt.xlabel('y position')
                plt.ylabel('y force')
                # print(y[i], F_g(data)[1])
            elif y_range[0] == y_range[1]:
                plt.plot(x[i], radiation_force_1(data)[0], 'g.')
                # print(x[i], F_s())
                plt.title('x position respect to force')
                plt.xlabel('x position')
                plt.ylabel('x force')
        plt.grid()
        plt.show()

    def inter2(x_range, y_range): # plot the force (calculated by newly developed method) respect to position
        x = np.linspace(x_range[0], x_range[1], 300)
        y = np.linspace(y_range[0], y_range[1], 300)
        plt.figure()
        for i in range(len(x)):
            input.droplet_pos_e = [x[i], y[i], 0, 0, 0, 0]
            # data, incident_angle, tilt = filter()
            print(i)
            rays, incident_angle, tilt = filter()
            if x_range[0] == x_range[1]:
                plt.plot(y[i], radiation_force_2(rays, incident_angle, tilt)[1], 'k.')
                plt.title('y position respect to force')
                plt.xlabel('y position')
                plt.ylabel('y force')
                # print(y[i], F_g(data)[1])
            elif y_range[0] == y_range[1]:
                plt.plot(x[i], radiation_force_2(rays, incident_angle, tilt)[0], 'g.')
                # print(x[i], F_s())
                plt.title('x position respect to force')
                plt.xlabel('x position')
                plt.ylabel('x force')
        plt.grid()
        plt.show()

    def inter_inter(x_range, y_range): # plot comparison between the two methods
        x = np.linspace(x_range[0], x_range[1], 300)
        y = np.linspace(y_range[0], y_range[1], 300)
        plt.figure()
        force1 = []
        force2 = []

        if x_range[0] == x_range[1]:
            for i in range(len(x)):
                print(i)
                rays, incident_angle, tilt = filter()
                input.droplet_pos_e = [x[i], y[i], 0, 0, 0, 0]
                force1.append(radiation_force_1(rays)[1])
                force2.append(radiation_force_2(rays, incident_angle, tilt)[1])

            plt.plot(y, force1, 'k.', label='Force calculated by method 1')
            plt.plot(y, force2, 'r.', label='Force calculated by method 2')
            plt.title('y position respect to force')
            plt.xlabel('y position(um)')
            plt.ylabel('y force(N)')

        elif y_range[0] == y_range[1]:
            for i in range(len(x)):
                input.droplet_pos_e = [x[i], y[i], 0, 0, 0, 0]
                print(i)
                rays, incident_angle, tilt = filter()
                force1.append(radiation_force_1(rays)[0])
                force2.append(radiation_force_2(rays, incident_angle, tilt)[0])
            plt.plot(x, force1, 'k.', label='Ashkin method')
            plt.plot(x, force2, 'r.', label='New code method')
            # plt.title('x position respect to force')
            plt.xlabel('x position(um)')
            plt.ylabel('x force(N)')
        plt.legend()
        plt.grid()
        plt.show()

    def inter_ellipse_shape(x_range, y_range):
        # plot the force (calculated by newly developed method) respect to position with different shape of ellipse
        x = np.linspace(x_range[0], x_range[1], 300)
        y = np.linspace(y_range[0], y_range[1], 300)
        plt.figure()

        major_axis = [10E-6,20E-6,30E-6,40E-6]
        colour = ['r','k','b','g']
        if x_range[0] == x_range[1]:
            for k in range(len(major_axis)):
                force = []
                input.a=major_axis[k]
                for i in range(len(x)):
                    print(i)
                    rays, incident_angle, tilt = filter()
                    input.droplet_pos_e = [x[i], y[i], 0, 0, 0, 0]
                    force.append(radiation_force_2(rays, incident_angle, tilt)[1])

                plt.plot(y, force, colour[k]+'.', label='major axis={} um'.format(int(major_axis[k]*1E6)))
            plt.title('y position respect to force')
            plt.xlabel('y position(um)')
            plt.ylabel('y force(N)')

        elif y_range[0] == y_range[1]:
            for k in range(len(major_axis)):
                force = []
                input.a=major_axis[k]
                for i in range(len(x)):
                    print(i)
                    rays, incident_angle, tilt = filter()
                    input.droplet_pos_e = [x[i], y[i], 0, 0, 0, 0]
                    force.append(radiation_force_2(rays, incident_angle, tilt)[0])

                plt.plot(x, force, colour[k]+'.', label='major axis={} um'.format(int(major_axis[k]*1E6)))
            plt.title('x position respect to force')
            plt.xlabel('x position(um)')
            plt.ylabel('x force(N)')

        plt.legend(loc='upper right')
        plt.grid()
        plt.show()

    def inter_ellipse_rotate(x_range, y_range):
        # plot the force (calculated by newly developed method) respect to rotational angle
        x = np.linspace(x_range[0], x_range[1], 300)
        y = np.linspace(y_range[0], y_range[1], 300)
        plt.figure()

        rotation = np.linspace(0, np.pi/2, 5)
        colour = ['r','k','b','g', 'y']
        if x_range[0] == x_range[1]:
            for k in range(len(rotation)):
                force = []
                # input.droplet_pos_e[2] = rotation[k]
                for i in range(len(x)):
                    rays, incident_angle, tilt = filter()
                    print(i)
                    input.droplet_pos_e = [x[i], y[i], rotation[k], 0, 0, 0]
                    force.append(radiation_force_2(rays, incident_angle, tilt)[1])

                plt.plot(y, force, colour[k]+'.', label='rotation={} degree'.format(int(rotation[k]*180/np.pi)))
            plt.title('y position respect to force')
            plt.xlabel('y position(um)')
            plt.ylabel('y force(N)')

        elif y_range[0] == y_range[1]:
            for k in range(len(rotation)):
                force = []
                # input.droplet_pos_e[2] = rotation[k]
                for i in range(len(x)):
                    rays, incident_angle, tilt = filter()
                    print(i)
                    input.droplet_pos_e = [x[i], y[i], rotation[k], 0, 0, 0]
                    force.append(radiation_force_2(rays, incident_angle, tilt)[0])
                plt.plot(x, force, colour[k]+'.', label='rotation={} degree'.format(int(rotation[k]*180/np.pi)))

            plt.title('x position respect to force')
            plt.xlabel('x position(um)')
            plt.ylabel('x force(N)')

        plt.legend(loc='upper right')
        plt.grid()
        plt.show()

    def inter_ellipse_rotate2(rot_range):
        # plot the torque respect to rotational angle
        angle = np.linspace(rot_range[0], rot_range[1], 300)
        plt.figure()
        force = []
        for i in range(len(angle)):
            rays, incident_angle, tilt = filter()
            print(i)
            input.droplet_pos_e[2] = angle[i]
            # print(input.droplet_pos_e)
            # force.append(radiation_force_2()[1])
            force.append(Torque.total_torque(rays, incident_angle, tilt, 1))
            print(Torque.total_torque(rays, incident_angle, tilt, 1))

        plt.plot(angle, force, 'k.')
        plt.title('torque vs rotation at [1000,0]')
        plt.xlabel('angle')
        plt.ylabel('torque')
        # plt.legend(loc='upper right')
        plt.grid()
        plt.show()
        print(force)



    start = timeit.default_timer()

    # inter1([900E-6, 900E-6], [-40E-6, 40E-6])
    # inter1([800E-6, 1000E-6], [0,0])

    # inter2([900E-6, 900E-6], [-40E-6,40E-6])
    # inter2([270E-6, 1000E-6], [0, 0])

    # inter_inter([900E-6, 900E-6], [-40E-6,40E-6])
    # inter_inter([270E-6, 1200E-6], [0,0])

    # inter_ellipse_shape([900E-6, 900E-6], [-40E-6, 40E-6])
    # inter_ellipse_shape([270E-6, 1200E-6], [0,0])

    # inter_ellipse_rotate([900E-6, 900E-6], [-40E-6, 40E-6])
    # inter_ellipse_rotate([270E-6, 1200E-6], [0,0])

    # inter_ellipse_rotate2([0, 2*np.pi])

    rays, incident_angle, tilt = filter()
    radiation_force_2(rays, incident_angle, tilt)

    stop = timeit.default_timer()
    print('Time: ', stop - start)

