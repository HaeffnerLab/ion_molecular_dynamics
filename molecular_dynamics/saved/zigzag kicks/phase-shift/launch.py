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

starting_positions[:, 2] = equil.get_positions(p.number_ions, p.f_z, p)

# starting_x,starting_y,starting_z = np.load('21_300000.0000.npy')#totally linear
# starting_x,starting_y,starting_z = np.load('21_310000.0000.npy')#still linear
# starting_x,starting_y,starting_z = np.load('21_318000.0000.npy')#still linear but borderline
# starting_x,starting_y,starting_z = np.load('21_319000.0000.npy')#smallest zigzag  
# starting_x,starting_y,starting_z = np.load('21_330000.0000.npy')#zigzag 
# starting_x,starting_y,starting_z = np.load('21_350000.0000.npy')#totally zigzag
# starting_x,starting_y,starting_z = np.load('21_400000.0000.npy')#deep zigzag 

# pyplot.figure()
# pyplot.plot(starting_z, starting_x, 'o')
# pyplot.ylim(pyplot.xlim())
# pyplot.show()
amplitude = 50e-9
# starting_positions[:, 0] = starting_x
# starting_positions[:, 1] = p.number_ions * [0]
# starting_positions[:, 2] = starting_z
starting_positions[7, 0] = amplitude

print starting_positions

positions,excitations = simulation.simulation(starting_positions, starting_velocities, random_seeding = 0)

larger_amplitude = 10**1 * amplitude
starting_positions[7,0] = larger_amplitude
positions_larger,excitations_larger = simulation.simulation(starting_positions, starting_velocities, random_seeding = 0)


# velocities = np.diff(positions, axis = 1) / p.timestep
# print 'vel'
time_axis = np.arange(positions.shape[1]) * p.timestep * 10**6
# print 't'
# energy_left = velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
# print 'eleft'
# energy_right = velocities[-1, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
# print 'eright'
# energy_center = velocities[10, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
# print 'center'
# 
# chunksize = 3000
# numchunks = energy_left.size // chunksize 
# energy_left = energy_left[:chunksize*numchunks].reshape((-1, chunksize))
# energy_right = energy_right[:chunksize*numchunks].reshape((-1, chunksize))
# energy_center = energy_center[:chunksize*numchunks].reshape((-1, chunksize))
# time_axis_binned = time_axis[:chunksize*numchunks].reshape((-1, chunksize))
# energy_left = energy_left.mean(axis=1)
# energy_right = energy_right.mean(axis=1)
# energy_center = energy_center.mean(axis=1)
# time_axis_binned = time_axis_binned.mean(axis=1)


pyplot.figure()
pyplot.plot(time_axis, positions[7, :, 0] / amplitude, label = 'left ion')
pyplot.plot(time_axis, 2 * positions_larger[7, :, 0] / larger_amplitude, label = 'left ion large displacement')
# pyplot.plot(time_axis, 10**6 * positions[-1, :, 0], label = 'right ion')
pyplot.title(r"5 ions, x - kick propagation, f_z = {} KHz".format(p.f_z / 10.**3))
pyplot.xlabel(r'$\mu s$')
pyplot.ylabel(r'Displacement along x, $\mu m$')
pyplot.tick_params(axis='both')
pyplot.legend()

# pyplot.figure()
# offset = 1
# max_energy = energy_center.max()
# energy_left = energy_left / max_energy
# energy_right = energy_right / max_energy + offset
# energy_center = energy_center / max_energy + 2 * offset
# pyplot.plot(time_axis_binned, energy_left, 'b', label = 'left ion')
# pyplot.plot(time_axis_binned, energy_right, 'r', label = 'right ion')
# pyplot.plot(time_axis_binned, energy_center, 'k', label = 'center ion')
# pyplot.legend()

# pyplot.figure()
# freq = np.fft.fftfreq(time_axis_binned.size, chunksize*p.timestep)
# fft = np.fft.fft(energy_left)
# total_fft = np.abs(fft)
# fft = np.fft.fft(energy_right)
# total_fft += np.abs(fft)
# #not plotting the DC term, only plotting f>0
# 
# pos_freqs = np.where(freq > 0)
# pyplot.plot(freq[pos_freqs] /10**3, total_fft[pos_freqs])
# pyplot.xlabel('Frequency KHz')
# pyplot.title('FFT', fontsize = 20)


pyplot.show()