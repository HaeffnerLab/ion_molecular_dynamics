import labrad.units as U
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

cdef void calculate_acceleration(double [:, :] position, double [:, :] current_acceleration, double time):
    '''
    given the current position, computes the current acceleration and fills in the current_acceleration array
    '''
    cdef int i = 0
    for i in range(NUMBER_IONS):
        current_acceleration[i, 0] = (-1/2.0) *( (sqrt(4) * W_DRIVE * sqrt(W_X ** 2 + W_Y **2 + W_Z**2)) * position[i, 0] * cos (W_DRIVE * time) + (- W_Z**2 + W_X**2 - W_Y**2) *  position[i, 0])  
        current_acceleration[i, 1] = (-1/2.0) *( - (sqrt(4) * W_DRIVE *  sqrt(W_X ** 2 +W_Y **2 + W_Z**2)) * position[i, 1] * cos (W_DRIVE * time) + (- W_Z**2 - W_X**2 + W_Y**2) *  position[i, 1]) 
        current_acceleration[i, 2] = (-1/2.0) *( 2 * W_Z**2 *  position[i, 2])

cdef void do_verlet_integration(double [:, :, :] positions):
    cdef double[:, :] current_position = positions[:,0,:]
    initial_acceleration =  np.zeros((NUMBER_IONS, 3))
    cdef double[:, :] current_acceleration = initial_acceleration
    cdef int i
    cdef int j
    cdef int k
    cdef double current_time
    cdef double next_print_progress = 1
    for i in range(2, TOTAL_STEPS):
        #print progress update
        if next_print_progress / 100.0 < float(i) / TOTAL_STEPS:
            print 'PROGRESS: {} %'.format(next_print_progress)
            next_print_progress = next_print_progress + 1
        current_time = i * TIMESTEP
        calculate_acceleration(current_position, current_acceleration, current_time)
        #cycle over ions
        for j in range(NUMBER_IONS):
            #cycle over coordinates
            for k in range(3):
                positions[j, i, k] = 2 * positions[j , i - 1, k] -  positions[j, i - 2, k] + current_acceleration[j, k] * TIMESTEP**2
        current_position = positions[:,i,:]
    
def simulation(starting_position):
    assert starting_position.shape == (NUMBER_IONS, 3), "Incorrect starting position format"
    positions = np.zeros((NUMBER_IONS, TOTAL_STEPS, 3))
    positions[:, 0, :] = starting_position
    positions[:, 1, :] = starting_position
    cdef double [:, :, :] positions_view = positions
    do_verlet_integration(positions_view)
    return positions