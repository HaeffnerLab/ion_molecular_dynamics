'''
python setup.py build_ext --inplace
'''

from simulation import simulator
import numpy as np
from simulation_parameters import simulation_parameters
from equilbrium_positions import equilibrium_positions as equil
import time


p = simulation_parameters()
simulation = simulator(p)

starting_positions = np.zeros((p.number_ions, 3))
starting_velocities = np.zeros((p.number_ions, 3))
starting_positions[:, 0] = p.number_ions * [0]
starting_positions[:, 1] = p.number_ions * [0]
starting_positions[:, 2] = equil.get_positions(p.number_ions, p.f_z, p)
starting_positions[0, 0] = 1e-6
start_time = t = time.time()
output,excitations = simulation.simulation(starting_positions, starting_velocities, random_seeding = 0)
finish_time = time.time()
print 'Simulation took {0:.2f} seconds'.format(finish_time - start_time)
# print output[0,:]

from matplotlib import pyplot
time_axis = np.arange(output.shape[1]) * p.timestep * 10**6
# for i in range(p.number_ions):
#     pyplot.plot(time_axis, output[i, :, 0])
pyplot.plot(time_axis, 10**6 * output[0, :, 0], label = 'left ion')
pyplot.plot(time_axis, 10**6 * output[-1, :, 0], label = 'right ion')
pyplot.title(r"5 ions, x - kick propagation, f_z = {} KHz".format(p.f_z / 10.**3), fontsize = 30)
pyplot.xlabel(r'$\mu s$', fontsize = 30)
pyplot.ylabel(r'Displacement along x, $\mu m$', fontsize = 30)
pyplot.tick_params(axis='both', labelsize = 20)
pyplot.xlim(0,200)
pyplot.legend(fontsize = 20)

pyplot.show()