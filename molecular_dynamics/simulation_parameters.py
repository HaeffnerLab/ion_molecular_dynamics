from __future__ import division

class simulation_parameters(object):
    
    number_ions = 1
    ion_amu_mass = 40
    #trap frequencies
    f_drve = 30.0 * 10**6#Hz
    f_x = 4.0 * 10**6#MHz
    f_y = 3.0 * 10**6#Hz
    f_z = 0.1 * 10**6#Hz
    
    simulation_duration = 0.0001#seconds
    timestep = (1 / f_drve) /100#seconds
    total_steps = int(simulation_duration / timestep)