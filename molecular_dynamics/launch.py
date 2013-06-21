'''
python setup.py build_ext --inplace
'''

import simulation
import numpy as np
from simulation_parameters import simulation_parameters as p
from equilbrium_positions import equilibrium_positions as equil
import time

starting_positions = np.zeros((p.number_ions, 3))

starting_positions[:, 0] = p.number_ions * [0]
starting_positions[:, 1] = p.number_ions * [0]
starting_positions[:, 2] = equil.get_positions(p.number_ions, p.f_z)
starting_positions[0, 0] = 10**-5
start_time = t = time.time()
output = simulation.simulation(starting_positions, random_seeding = 0)
finish_time = time.time()
print 'Simulation took {0:.2f} seconds'.format(finish_time - start_time)
# print output[0,:]

from matplotlib import pyplot
time_axis = np.arange(output.shape[1]) * p.timestep * 10**6
for i in range(p.number_ions):
    pyplot.plot(time_axis, output[i, :, 0])

pyplot.show()