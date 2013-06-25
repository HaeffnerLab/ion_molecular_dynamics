from __future__ import division
import numpy as np

class simulation_parameters(object):
    
    number_ions = 3
    #trap frequencies
    f_drve = 30.0 * 10**6#Hz
    f_x = 4.0 * 10**6#MHz
    f_y = 3.0 * 10**6#Hz
    f_z = 0.2 * 10**6#Hz
    #simulation parameters
    damping = 0 #optional velocity damping, useful for finding equlibrium positions
    simulation_duration = 0.002#seconds
    timestep = (1 / f_drve) /100#seconds
    total_steps = int(simulation_duration / timestep)
    #ion parameters
    mass = 40 * 1.6605402e-27 #40 amu in kg
    coulomb_coeff = 2.30707955552e-28 # k =  U.e**2 / (4.0 * U.pi * U.eps0) 
    hbar = 1.05457266913e-34
    transition_gamma = (1 / (7.1 * 10**-9)) #Gamma = 1 / Tau
    #laser
    saturation =1.0
    laser_detuning = -.5 * transition_gamma
    laser_direction = np.array([1,1,1]); laser_direction = laser_direction / np.sqrt(np.sum(laser_direction))#normalized
    transition_k_mag =  2 * np.pi / (396.847 * 10**-9) 