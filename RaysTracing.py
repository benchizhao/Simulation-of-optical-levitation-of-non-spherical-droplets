# -*- coding: utf-8 -*-
"""
Time: 2020/Jul/16

@author: Benchi Zhao
Python version: 3.7
Functionality: This module trace all individual rays and plot the diagram. Also plot the intensity of cross-section.
Dependence: Need import 'input', all parameters of the experiment is saved in that file 'input'.
"""

from RayTracing import Ray
import input
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import timeit

def gaussian(y, sigma):
    '''
    :param y: the y position describes the distance between beam axis and the the ray
    :return: the gaussian value at y
    '''
    return 1 / (np.sqrt(2 * np.pi)*sigma) * np.exp(-y ** 2 / (2*sigma**2))

def intensity(position):
    '''
    Generate the intensity of bundle of rays
    :param position: list
        there are input.no_of_rays values in the list, give the intensity of each rays at the beginning
    :return: list
        normalised intensity list
    '''
    I = []
    for i in range(len(position)):
        intens = gaussian(position[i], input.sigma)
        I.append(intens)
    normalise_I = I / sum(I)
    return normalise_I

def bundle():
    '''
    Trace a bundle of rays, the number of rays i determined by 'input' file, no_of_rays.
    the function returns a list called rays. The shape of the list is (no_of_rays, x, 4)
     --x is 8 if the ray do not interact with droplet;
     --x is 10 is the ray interact with droplet.
     :return: list
        rays: list of ray_state, the shape is (no_of_rays, x, 4) the '4' represent [x,y,theta,intensity]
        incident_angle: list of the incident angle between ray and curved surface(lens/ droplet), the shape is (no_of_rays, x)
                        if the rays interact with droplet, x is 2
                        if the rays do not interact with droplet, x is 4
        tilt: list of the tilt angle(tangent) of the curved surface(lens/ droplet), the shape is (no_of_rays, x)
              if the rays interact with droplet, x is 2
              if the rays do not interact with droplet, x is 4
              tilt angle is 0 if the tangent of the curvature is vertical;
              tilt angle is positive if the slope of the tangent is negative;
              tilt angle is negative if the slope of the tangent is positive.
    '''
    rays = []
    incident_angle = []
    tilt = []
    list_of_rays = np.linspace(-input.width, input.width, input.no_of_rays)

    if input.droplet_shape in ['circle', 'c']:
        for i in range(input.no_of_rays):
            r = Ray()
            r.trigger(list_of_rays[i], intensity(list_of_rays)[i])
            r.lens(f=input.lens_f, position=input.lens_pos, thick=input.len_thickness)
            r.circle_droplet(input.radius)
            rays.append(r.ray_state)
            incident_angle.append(r.incident_angle)
            tilt.append(r.tilt)

    elif input.droplet_shape in ['ellipse', 'e']:
        for i in range(input.no_of_rays):
            r = Ray()
            r.trigger(list_of_rays[i], intensity(list_of_rays)[i])
            r.lens(f=input.lens_f, position=input.lens_pos, thick=input.len_thickness)
            r.ellipse_droplet(a=input.a, b=input.b, center_point=input.droplet_pos_e[:2], theta=input.droplet_pos_e[2])
            rays.append(r.ray_state)
            incident_angle.append(r.incident_angle)
            tilt.append(r.tilt)
    return rays, incident_angle, tilt

def make_table(data):
    '''
       make_table(data)
           Show the ray states in table form,  making it clean.

       Parameters
       ----------
       data: list
           all of the ray state (return of function bundle).

       Returns
       -------
           Show the data in table form.

    '''
    for i in range(input.no_of_rays):
        data_1 = data[i]
        df_1 = pd.DataFrame(data_1, columns=['x', 'y', 'theta', 'intensity'])
        print('Ray', i + 1)
        print(df_1)

def plot_bundle_c(data):
    '''
    PLot the bundle of rays passing through the lens and circular droplet.
    :param data:
        The ray state of all rays
    '''
    plt.figure('Ray Tracing')

    for i in range(input.no_of_rays):
        for j in range(len(data[i])-1):
            x = np.linspace(data[i][j][0], data[i][j+1][0])
            y = np.linspace(data[i][j][1], data[i][j+1][1])
            plt.plot(x, y, 'k-', linewidth=data[i][j][3]*20)

        # plt.xlim(0, max_x)

    theta = np.linspace(0, 2 * np.pi)
    x_c = input.radius * np.cos(theta) + input.droplet_pos_e[0]
    y_c = input.radius * np.sin(theta) + input.droplet_pos_e[1]
    plt.plot(x_c, y_c)

    # # plot lens
    theta1 = np.linspace(3 * np.pi / 4, 5 * np.pi / 4)
    theta2 = np.linspace(- np.pi / 4, np.pi / 4)
    x_1 = input.lens_f * np.cos(theta1) + input.lens_pos + input.lens_f - input.len_thickness / 2
    x_2 = input.lens_f * np.cos(theta2) + input.lens_pos - input.lens_f + input.len_thickness / 2
    y_1 = input.lens_f * np.sin(theta1)
    y_2 = input.lens_f * np.sin(theta2)
    plt.plot(x_1, y_1)
    plt.plot(x_2, y_2)
    # plt.grid()
    plt.ylim(-50E-6, 50E-6)
    plt.show()

def plot_bundle_e(data):
    '''
       PLot the bundle of rays passing through the lens and circular droplet.
    :param data: float
        The ray state of all rays
    '''
    plt.figure('Ray Tracing')

    for i in range(input.no_of_rays):
        for j in range(len(data[i]) - 1):
            x = np.linspace(data[i][j][0], data[i][j + 1][0])
            y = np.linspace(data[i][j][1], data[i][j + 1][1])
            plt.plot(x, y, 'k-', linewidth=data[i][j][3]*20)#
    # plot droplet
    theta = np.linspace(0, 2 * np.pi, 100)
    a = input.a
    b = input.b
    e_rot_angle = input.droplet_pos_e[2]
    center_point = input.droplet_pos_e[:2]
    x = a * np.cos(theta)
    y = b * np.sin(theta)
    x_e = np.cos(e_rot_angle) * x - np.sin(e_rot_angle) * y + center_point[0]
    y_e = np.sin(e_rot_angle) * x + np.cos(e_rot_angle) * y + center_point[1]
    plt.plot(x_e, y_e)

    # # plot lens
    theta1 = np.linspace(3 * np.pi / 4, 5 * np.pi / 4)
    theta2 = np.linspace(- np.pi / 4, np.pi / 4)
    x_1 = input.lens_f * np.cos(theta1) + input.lens_pos + input.lens_f - input.len_thickness / 2
    x_2 = input.lens_f * np.cos(theta2) + input.lens_pos - input.lens_f + input.len_thickness / 2
    y_1 = input.lens_f * np.sin(theta1)
    y_2 = input.lens_f * np.sin(theta2)
    plt.plot(x_1, y_1, color='grey')
    plt.plot(x_2, y_2, color='grey')
    # plt.grid()
    plt.ylim(-50E-6, 50E-6)
    plt.show()

def detector(x_position, data):
    '''
    Detect the intensity of rays at position x.
    :param x_position: float
        The position we would like to detect the bundle of rays
    :param data: list
        The list contains all ray state data
    :return y_values: list
        positions of rays in y-axis.
    :return intens_value: list
        the corresponding intensity values at position y.

    '''
    y_values = []
    intns_values = []
    for i in range(input.no_of_rays):
        refract_x = np.array(data[i])[:, 0]
        refract_y = np.array(data[i])[:, 1]
        refract_theta = np.array(data[i])[:, 2]
        refract_intensity = np.array(data[i])[:, 3]
        for j in range(len(refract_x)-1):
            if x_position > refract_x[j] and x_position < refract_x[j + 1]:
                y = np.tan(refract_theta[j]) * (x_position - refract_x[j]) + refract_y[j]
                y_values.append(y)
                intns_values.append(refract_intensity[j])
    return y_values, intns_values

def plot_detector(x_positions, data):
    '''
    :param x_positions: list
        input all x-positions where the intensity that we are interested at.
    :param data: list
        The data is the return(result) of function bundle()
    :return: output the figure of intensity distribution of the beam
    '''
    all_pos = []
    all_ints = []
    for k in range(len(x_positions)):
        y_values, intns_values = detector(x_positions[k], data)
        num = round((max(y_values) - min(y_values)) /2 * 1E6)
        n, bins, patches = plt.hist(y_values, int(num))
        index = 0
        binned = []
        b = [int(i) for i in n]
        for i in range(0, len(b), 1):
            c = intns_values[index:index + b[i]]
            index = index + b[i]
            binned.append(c)
        pos = []
        for i in range(len(bins)-1):
            ave = (bins[i] + bins[i+1])/2
            pos.append(ave)
        ints = [sum(i) for i in binned]
        all_pos.append(pos)
        all_ints.append(ints)

    plt.figure('Intensity distribution of the beam')
    for i in range(len(x_positions)):
        plt.plot(all_pos[i], all_ints[i], label="x = %f (um)" %(x_positions[i]*10E6))

    plt.legend(fontsize=12)
    plt.xlabel('Y-position')
    plt.ylabel('Intensity')
    plt.title('Intensity distribution of the beam')
    plt.show()

if __name__ == '__main__':
    '''Trace all individual rays'''
    start = timeit.default_timer()

    data, incident_angle, tilt = bundle() # achieve the ray state
    make_table(data=data) # show the ray state in table form

    if input.droplet_shape in ['circle', 'c']: # plot the rays if the target is circle
        plot_bundle_c(data)
    elif input.droplet_shape in ['ellipse', 'e']: # plot the rays if the target is ellipse
        plot_bundle_e(data)

    end = timeit.default_timer()
    print('time:', end - start)


    plot_detector([500E-6, 700E-6], data) # two plots will showed here, the first one is the intensitu distribution, and please ignore the second.