from simulation_parameters import simulation_parameters as p
from libc.math cimport cos, sqrt
import numpy as np

cdef int TOTAL_STEPS = p.total_steps
cdef int NUMBER_IONS = p.number_ions
cdef double TIMESTEP = p.timestep
cdef double W_DRIVE = p.f_drve * 2 * np.pi
cdef double W_X = p.f_x * 2 * np.pi
cdef double W_Y = p.f_y * 2 * np.pi
cdef double W_Z = p.f_z * 2 * np.pi
cdef double MASS = p.mass
cdef double COULOMB_K = p.coulomb_k
cdef double VEL_DAMPING = p.damping


cdef void calculate_acceleration(double [:, :] position, double [:, :] velocity, double [:, :] current_acceleration, double time, double random_float):
    '''
    given the current position, computes the current acceleration and fills in the current_acceleration array
    '''
    cdef int i = 0
    cdef int j = 0
    cdef double dx, dy, dz, distance_sq
    cdef double Fx, Fy, Fz
    #acceleration due to the trap
    for i in range(NUMBER_IONS):
        current_acceleration[i, 0] = ( (1/2.)*(-W_X**2 + W_Y**2 + W_Z**2) - W_DRIVE * sqrt(W_X**2 + W_Y**2 + W_Z**2) * cos (W_DRIVE * time)) * position[i, 0]
        current_acceleration[i, 1] = ( (1/2.)*( W_X**2 - W_Y**2 + W_Z**2) + W_DRIVE * sqrt(W_X**2 + W_Y**2 + W_Z**2) * cos (W_DRIVE * time)) * position[i, 1]
        current_acceleration[i, 2] = - W_Z**2 *  position[i, 2]
        #acceleration due to the coulombic repulsion
    for i in range(NUMBER_IONS):
        for j in range(i + 1, NUMBER_IONS):
            dx = position[i, 0] - position[j, 0]
            dy = position[i, 1] - position[j, 1]
            dz = position[i, 2] - position[j, 2]
            distance_sq = dx**2 + dy**2 + dz**2
            if distance_sq == 0:
                raise Exception("Distance between ions is 0")
            Fx = COULOMB_K * dx / (distance_sq)**(3./2.)
            Fy = COULOMB_K * dy / (distance_sq)**(3./2.)
            Fz = COULOMB_K * dz / (distance_sq)**(3./2.)
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

cdef void do_verlet_integration(double [:, :, :] positions, double [:] random_floats):
    cdef double[:, :] current_position = positions[:,0,:]
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
        calculate_acceleration(current_position, current_velocity, current_acceleration, current_time, random_floats[i])
        #cycle over ions
        for j in range(NUMBER_IONS):
            #cycle over coordinates
            for k in range(3):
                positions[j, i, k] = 2 * positions[j , i - 1, k] -  positions[j, i - 2, k] + current_acceleration[j, k] * TIMESTEP**2
                current_velocity[j, k] = (positions[j, i, k] - positions[j, i - 2, k]) / (2 * TIMESTEP)
        current_position = positions[:,i,:]
    
def simulation(starting_position, random_seeding = None):
    assert starting_position.shape == (NUMBER_IONS, 3), "Incorrect starting position format"
    if random_seeding is not None:
        np.random.seed(random_seeding)
    #precalcualte random floats we will be using
    random_floats = np.random.random(TOTAL_STEPS)
    positions = np.zeros((NUMBER_IONS, TOTAL_STEPS, 3))
    positions[:, 0, :] = starting_position
    positions[:, 1, :] = starting_position
    cdef double [:, :, :] positions_view = positions
    cdef double [:] random_float_view = random_floats
    do_verlet_integration(positions_view, random_float_view)
    return positions