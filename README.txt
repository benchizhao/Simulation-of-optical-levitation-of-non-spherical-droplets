SIMULATION OF OPTICAL LEVIATION OF NON-SPHERICAL DROPLETS.
=========================================================
#####################
AUTHOR: 
-------
	BENCHI ZHAO
	bz819@ic.ac.uk
#####################

What does the code do?
------------------------
This code models the movement of spherical and non-spherical droplets inside a Gaussian beam 
laser. The code calculates the optical forces on the target and updates the position and 
velocity by solving Newton’s 2nd law.


Dependencies and Installation:
------------------------------
The code uses and has been tested with the following packages:
Python Version: 3.7 
Numpy Version: 1.15.4
Scipy Version: 1.1.0
Matplotlib Version: 3.0.2
Pandas Version: 0.23.4
Math Version: 1.1.0

These can be installed using a Python distribution such as Anaconda:
https://www.anaconda.com/
Required packages installed by selecting the Python 3.7 version for your platform.
The code has been run extensively on Windows. The code has not been run on IOS or Linux.



###########################################
There are 12 .py documents in this project.
###########################################

# input.py
----------
This module determines the parameters used in the simulation, e.g. laser power, droplet 
size, density or air, etc.

# ellipse.py
------------
This module is to finding the interaction between straight line and ellipse, the lintersection
points between dtraight line and ellpse; find the slope of the line that tangent to the ellipse.


# RayTracing.py
---------------
This module is to tracing an individual ray passing the lens and droplet, and plot the 
path of the ray.

# RaysTracing.py
---------------
This module is to tracing all individual rays passing the lens and droplet, and plot the 
path of the ray. And could plot the intensity at different position.

# Radiation_Force.py
--------------------
This module calculates the raidation force induced by laser by using both Askin's method 
and newly developed method.

# Drag_force.py
---------------
This module calculates the drag force (air resistance) when droplet moves in the air.

# Buoyancy.py
-------------
This module calculates the buoyancy force in the system, which is a constant force
encounters the gravity.

# Torque.py
------------
This module is to calculating the momentum of inertia off ellipse and torque caused by 
radiation pressure and air resistance.

# movement.py
-------------
This module is solving the trajectories of the droplet by ODE solvers.

# GPU programming.py
--------------------
This module is a code to test how fast the GPU programming could fast the code.

# Multiprocesses.py
-----------------
This module is a code to test how fast the Mulitiprocess could fast the code.

# showAshkin.py
---------------
This module is repeating the results achieved by Ahskin.

####################
How to run the code?
####################


1. Set up the parameters of the droplet in the input.py file. For example, the shape of droplet(circle/ellipse), the radius/major-axis, 
the power of the laser, the initial poition of droplets, except the power mode, which is set in movement.py.

2. There are three modules can be ignored: GPU programming.py; show Ashkin.py; 
Mulitiprocesses.py. These 3 modules are designed for test.

3. Test if other 9 modules works. The codes written under
if__name__=='__main__':
are test parts. They supposed to execute without bugs.
RayTracing.py generates figure of how one individual ray propagtes in the system.
RaysTracing.py generates figure of how a bundle of rays propagte in the system, and a beam intensity distribution figure.
Radiation_Force.py generates the figure of the relation between forces repects to the changing x or y.
Torque.py, Buoyancy.py and Drag_force.py generate the corresonding values.

4. After testing these codes, we should run the movement.py, which is the core part for solving the Newton's second law.
In this modules, forces and torques are prepared. In the function of diff_ellipse(d_list, t), we could change the laser power here,
the choices are constant power/noisy power, and the square function (the commented code), also, we can define any shape of the power
as we want (write your own code substitute the commemted code). The data generated in this module will be saved in test.txt file ,
the name can be changed. 
Then the figrue could plot how the centeral point of the droplet move on x and y axis respect to time. 










