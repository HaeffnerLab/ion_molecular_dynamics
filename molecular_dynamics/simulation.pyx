from simulation_parameters import simulation_parameters as p
from libc.math cimport cos, sin, sqrt, acos, M_PI
import numpy as np

cdef int TOTAL_STEPS = p.total_steps
cdef int NUMBER_IONS = p.number_ions
cdef double TIMESTEP = p.timestep
cdef double W_DRIVE = p.f_drve * 2 * np.pi
cdef double W_X = p.f_x * 2 * np.pi
cdef double W_Y = p.f_y * 2 * np.pi
cdef double W_Z = p.f_z * 2 * np.pi
cdef double MASS = p.mass
cdef double COULOMB_COEFF = p.coulomb_coeff
cdef double VEL_DAMPING = p.damping
cdef double HBAR = p.hbar
#laser parameters
cdef double GAMMA = p.transition_gamma
cdef double SATURATION = p.saturation
cdef double LASER_DETUNING = p.laser_detuning
cdef double [:] LASER_DIRECTION = p.laser_direction
cdef double TRANSITION_K_MAG = p.transition_k_mag

cdef void calculate_acceleration(double [:, :] position, double [:, :] velocity, double [:, :] current_acceleration, double time, double [:, :] random_floats, char[:] excitation):
    '''
    given the current position, computes the current acceleration and fills in the current_acceleration array
    '''
    cdef int i = 0
    cdef int j = 0
    cdef double dx, dy, dz, distance_sq
    cdef double Fx, Fy, Fz
    cdef double inst_detuning
    cdef double gamma_laser
    cdef double p_excited
    cdef double theta = 0 
    cdef double phi = 0
    #acceleration due to the trap
    for i in range(NUMBER_IONS):
        current_acceleration[i, 0] = ( (1/2.)*(-W_X**2 + W_Y**2 + W_Z**2) - W_DRIVE * sqrt(W_X**2 + W_Y**2 + W_Z**2) * cos (W_DRIVE * time)) * position[i, 0]
        current_acceleration[i, 1] = ( (1/2.)*( W_X**2 - W_Y**2 + W_Z**2) + W_DRIVE * sqrt(W_X**2 + W_Y**2 + W_Z**2) * cos (W_DRIVE * time)) * position[i, 1]
        current_acceleration[i, 2] = - W_Z**2 *  position[i, 2]
    #acceleration due to the coulombic repulsion
    for i in range(NUMBER_IONS):
        for j in range(i + 1, NUMBER_IONS):
            #the double for loop iterations over unique pairs
            dx = position[i, 0] - position[j, 0]
            dy = position[i, 1] - position[j, 1]
            dz = position[i, 2] - position[j, 2]
            distance_sq = dx**2 + dy**2 + dz**2
            if distance_sq == 0:
                raise Exception("Distance between ions is 0")
            Fx = COULOMB_COEFF * dx / (distance_sq)**(3./2.)
            Fy = COULOMB_COEFF * dy / (distance_sq)**(3./2.)
            Fz = COULOMB_COEFF * dz / (distance_sq)**(3./2.)
            #acceleartion is equal and opposite
            current_acceleration[i, 0] += Fx / MASS
            current_acceleration[i, 1] += Fy / MASS
            current_acceleration[i, 2] += Fz / MASS
            current_acceleration[j, 0] -= Fx / MASS
            current_acceleration[j, 1] -= Fy / MASS
            current_acceleration[j, 2] -= Fz / MASS
        #optional velocity damping
        current_acceleration[i, 0] += - VEL_DAMPING * velocity[i, 0]
        current_acceleration[i, 1] += - VEL_DAMPING * velocity[i, 1]
        current_acceleration[i, 2] += - VEL_DAMPING * velocity[i, 2]
    #acceleration due to laser interaction
    for i in range(NUMBER_IONS):
        inst_detuning = LASER_DETUNING - TRANSITION_K_MAG * ( LASER_DIRECTION[0] * velocity[i, 0] + LASER_DIRECTION[1] * velocity[i, 1]+ LASER_DIRECTION[2] * velocity[i, 2]) #Delta + k . v
        gamma_laser = SATURATION /  (1. + (2 * inst_detuning /  GAMMA)**2) * GAMMA / 2.
        if not excitation[i]:
            #if atom is not currently excited, calculate probability to get excited
            p_exc = gamma_laser * TIMESTEP
            if random_floats[i, 0] <= p_exc:
                #atom gets excited
                excitation[i] = 1
                #acceleration = (momentum change) / ( time * mass)
                current_acceleration[i, 0] += HBAR * TRANSITION_K_MAG * LASER_DIRECTION[0] / (TIMESTEP * MASS) 
                current_acceleration[i, 1] += HBAR * TRANSITION_K_MAG * LASER_DIRECTION[1] / (TIMESTEP * MASS) 
                current_acceleration[i, 2] += HBAR * TRANSITION_K_MAG * LASER_DIRECTION[2] / (TIMESTEP * MASS) 
            else:
                #atom stays in the ground state
                pass
        else:
            #if atom is currently excited, it can decay sponteneously or stimulated
            p_spon = GAMMA * TIMESTEP
            p_stim = gamma_laser * TIMESTEP
            if p_spon > 0.1 or p_stim > 0.1:
                raise Exception ("time step too small to deal with emission")
            if random_floats[i, 0] <= p_stim:
                #atom gets de-excited, stimualted
                excitation[i] = 0
                current_acceleration[i, 0] -= HBAR * TRANSITION_K_MAG * LASER_DIRECTION[0] / (TIMESTEP * MASS) 
                current_acceleration[i, 1] -= HBAR * TRANSITION_K_MAG * LASER_DIRECTION[1] / (TIMESTEP * MASS) 
                current_acceleration[i, 2] -= HBAR * TRANSITION_K_MAG * LASER_DIRECTION[2] / (TIMESTEP * MASS) 
            elif p_stim < random_floats[i, 0] <= (p_stim + p_spon):
                #atom gets de-excited, spontaneous, in a random direction
                excitation[i] = 0
                theta = 2 * M_PI * random_floats[i, 1]
                phi = acos(2 * random_floats[i, 2] - 1)
                current_acceleration[i, 0] += HBAR * TRANSITION_K_MAG * sin(theta) * cos(phi)   / (TIMESTEP * MASS) 
                current_acceleration[i, 1] += HBAR * TRANSITION_K_MAG * sin(theta) * sin(phi)   / (TIMESTEP * MASS) 
                current_acceleration[i, 2] += HBAR * TRANSITION_K_MAG * cos(theta)              / (TIMESTEP * MASS) 
            else:
                #atom stays excited
                pass
                

cdef void do_verlet_integration(double [:, :, :] positions, double [:, :, :] random_floats, char [:, :] excitations):
    cdef double[:, :] current_position = positions[:,0,:]
    initial_excitation = np.zeros(NUMBER_IONS, dtype = np.uint8)
    cdef char[:] current_excitation = initial_excitation
    initial_acceleration =  np.zeros((NUMBER_IONS, 3))
    cdef double[:, :] current_acceleration = initial_acceleration
    current_velocity_np = np.zeros((NUMBER_IONS, 3))
    cdef double[:, :] current_velocity = current_velocity_np
    cdef int i
    cdef int j
    cdef int k
    cdef double current_time
    cdef double next_print_progress = 0
    for i in range(2, TOTAL_STEPS):
        #print progress update
        if next_print_progress / 100.0 < float(i) / TOTAL_STEPS:
            print 'PROGRESS: {} %'.format(next_print_progress)
            next_print_progress = next_print_progress + 10
        current_time = i * TIMESTEP
        calculate_acceleration(current_position, current_velocity, current_acceleration, current_time, random_floats[i], current_excitation)
        #cycle over ions
        for j in range(NUMBER_IONS):
            #cycle over coordinates
            for k in range(3):
                positions[j, i, k] = 2 * positions[j , i - 1, k] -  positions[j, i - 2, k] + current_acceleration[j, k] * TIMESTEP**2
                current_velocity[j, k] = (positions[j, i, k] - positions[j, i - 2, k]) / (2 * TIMESTEP)
            excitations[i,j] = current_excitation[j]
        current_position = positions[:,i,:]
    
def simulation(starting_position, random_seeding = None):
    assert starting_position.shape == (NUMBER_IONS, 3), "Incorrect starting position format"
    if random_seeding is not None:
        np.random.seed(random_seeding)
    #precalcualte random floats we will be using
    random_floats = np.random.random((TOTAL_STEPS, NUMBER_IONS, 3))
    positions = np.zeros((NUMBER_IONS, TOTAL_STEPS, 3))
    excitations = np.zeros((TOTAL_STEPS, NUMBER_IONS), dtype = np.uint8)
    positions[:, 0, :] = starting_position
    positions[:, 1, :] = starting_position
    cdef double [:, :, :] positions_view = positions
    cdef double [:, :, :] random_float_view = random_floats
    cdef char [:,:] excitations_view = excitations
    do_verlet_integration(positions_view, random_float_view, excitations_view)
    return positions, excitations