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
starting_positions[0, 0] = 0
start_time = t = time.time()
positions, excitations = simulation.simulation(starting_positions, random_seeding = 15)
print positions.shape
velocities = np.diff(positions, axis = 1) / p.timestep
print velocities.shape


finish_time = time.time()
print 'Simulation took {0:.2f} seconds'.format(finish_time - start_time)
# print output[0,:]

from matplotlib import pyplot
time_axis = np.arange(positions.shape[1]) * p.timestep * 10**6
energy_x = velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
# for i in range(p.number_ions):
#     pyplot.plot(time_axis, output[i, :, 0])
# pyplot.plot(time_axis[:-1], velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x), label = 'left ion')
# pyplot.plot(time_axis, 10**6 * positions[-1, :, 0], label = 'right ion')
# pyplot.plot(time_axis, 10**6 * output[2, :, 0], label = 'center ion')
pyplot.title(r"5 ions, x - kick propagation, f_z = {} KHz".format(p.f_z / 10.**3), fontsize = 30)
pyplot.xlabel(r'$\mu s$', fontsize = 30)
pyplot.ylabel(r'Displacement along x, $\mu m$', fontsize = 30)
pyplot.tick_params(axis='both', labelsize = 20)
# pyplot.xlim(0,2000)
# pyplot.legend(fontsize = 20)

# pyplot.show()
print 'mean excitation', excitations.mean(axis = 0)

#downsampling the data to plot
chunksize = 1000
numchunks = energy_x.size // chunksize 
ychunks = energy_x[:chunksize*numchunks].reshape((-1, chunksize))
xchunks = time_axis[:chunksize*numchunks].reshape((-1, chunksize))
 
# Calculate the max, min, and means of chunksize-element chunks...
max_env = ychunks.max(axis=1)
min_env = ychunks.min(axis=1)
ycenters = ychunks.mean(axis=1)
xcenters = xchunks.mean(axis=1)
 
# Now plot the bounds and the mean...
pyplot.fill_between(xcenters, min_env, max_env, color='gray', 
                 edgecolor='none', alpha=0.5)
pyplot.plot(xcenters, ycenters)
pyplot.show()