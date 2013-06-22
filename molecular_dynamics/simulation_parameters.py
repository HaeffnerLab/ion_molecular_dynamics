from __future__ import division

class simulation_parameters(object):
    
    number_ions = 5
    #trap frequencies
    f_drve = 30.0 * 10**6#Hz
    f_x = 4.0 * 10**6#MHz
    f_y = 3.0 * 10**6#Hz
    f_z = 0.2 * 10**6#Hz
    #simulation parameters
    damping = 0 #optional velocity damping, useful for finding equlibrium positions
    simulation_duration = 0.001#seconds
    timestep = (1 / f_drve) /100#seconds
    total_steps = int(simulation_duration / timestep)
    #ion parameters
    mass = 40 * 1.6605402e-27 #40 amu in kg
    coulomb_coeff = 2.30707955552e-28 # k =  U.e**2 / (4.0 * U.pi * U.eps0) 