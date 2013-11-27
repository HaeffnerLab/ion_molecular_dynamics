'''
python setup.py build_ext --inplace
'''

from simulation import simulator
import numpy as np
from simulation_parameters import simulation_parameters
from equilbrium_positions import equilibrium_positions as equil
from matplotlib import pyplot

p = simulation_parameters()
simulation = simulator(p)

starting_positions = np.zeros((p.number_ions, 3))
starting_velocities = np.zeros((p.number_ions, 3))


starting_positions[:, 0] = 1e-6 * (np.random.rand(p.number_ions) - .5)
starting_positions[:, 1] = 1e-6 * (np.random.rand(p.number_ions) - .5)
starting_positions[:, 2] = (np.arange(p.number_ions) - 10) * 1e-6 

# equil.get_positions(p.number_ions, p.f_z, p)

positions,excitations = simulation.simulation(starting_positions, starting_velocities, random_seeding = 0)

final_x = positions[:,-1,0]
final_y = positions[:,-1,1]
final_z = positions[:,-1,2]
pyplot.figure()
pyplot.plot(positions[0,:,0]**2)
pyplot.yscale('log')

pyplot.figure()
pyplot.plot(final_z, final_x, 'o')
pyplot.ylim(pyplot.xlim())

order = np.argsort(final_z)
final_x = final_x[order]
final_y = final_y[order]
final_z = final_z[order]

np.save('{0:.0f}_{1:.4f}'.format(p.number_ions, p.f_z), (final_x,final_y,final_z))
pyplot.show()