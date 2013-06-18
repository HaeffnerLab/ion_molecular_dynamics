'''
python setup.py build_ext --inplace
'''

import simulation
import numpy as np
from simulation_parameters import simulation_parameters as p
import time

starting_positions = np.zeros((p.number_ions, 3))
starting_positions[:, 0] = np.array(0.000001)
# starting_positions[:, 2] = np.array([1])
start_time = t = time.time()
output = simulation.simulation(starting_positions)
finish_time = time.time()
print 'Simulation took {0:.2f} seconds'.format(finish_time - start_time)

print output.shape
# print output[0,:]

from matplotlib import pyplot
time_axis = np.arange(output.shape[1]) * p.timestep * 10**6
pyplot.plot(time_axis, output[0, :, 0], '-')
pyplot.show()