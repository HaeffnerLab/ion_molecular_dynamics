'''
python setup.py build_ext --inplace
'''
from simulation import simulator
import numpy as np
from simulation_parameters import simulation_parameters
from equilbrium_positions import equilibrium_positions as equil
from matplotlib import pyplot

doppler_average = np.zeros(599)
heating_average = np.zeros(599)
total_to_average = 10
for i in range(total_to_average):
    print 'ITERATION', i
    #starting at equilibrium positions
    p = simulation_parameters()
    starting_positions = np.zeros((p.number_ions, 3))
    starting_positions[:, 0] = p.number_ions * [0]
    starting_positions[:, 1] = p.number_ions * [0]
    starting_positions[:, 2] = equil.get_positions(p.number_ions, p.f_z, p)
    starting_velocities = np.zeros((p.number_ions, 3))
    #first do doppler cooling to randomize starting positions
    p.laser_detuning = -.5 * p.transition_gamma
    p.simulation_duration = 200e-6
    doppler_cooling_simulator = simulator(p)
    doppler_positions, doppler_excitations = doppler_cooling_simulator.simulation(starting_positions, starting_velocities, random_seeding = i)
    #now do laser heating
    p.laser_detuning = +0.5 * p.transition_gamma
    p.saturation = 1.0
    p.simulation_duration = 200e-6
    p.laser_direction = np.array([1.0, 0.0, 0.0])
    starting_positions = doppler_positions[:,-1,:]
    starting_velocities =  doppler_positions[:,-1,:] -  doppler_positions[:,-2,:]
    heating_simulator = simulator(p)
    heating_positions, heating_excitations = heating_simulator.simulation(starting_positions, starting_velocities, random_seeding = i)
    
    
    doppler_velocities = np.diff(doppler_positions, axis = 1) / p.timestep
    heat_velocities = np.diff(heating_positions, axis = 1) / p.timestep
    
    time_axis_doppler = np.arange(doppler_positions.shape[1]) * p.timestep * 10**6
    time_axis_heat = time_axis_doppler[-1] + np.arange(heating_positions.shape[1]) * p.timestep * 10**6
    doppler_energy_x = doppler_velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
    heat_energy_x = heat_velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
    # pyplot.title(r"5 ions, x - kick propagation, f_z = {} KHz".format(p.f_z / 10.**3), fontsize = 30)
    pyplot.xlabel(r'$\mu s$', fontsize = 30)
    pyplot.ylabel(r'Motional Quanta x', fontsize = 30)
    pyplot.tick_params(axis='both', labelsize = 20)
    # pyplot.legend(fontsize = 20)
     
    #downsampling the data to plot
    chunksize = 1000
    numchunks_doppler = doppler_energy_x.size // chunksize 
    ychunks_doppler = doppler_energy_x[:chunksize*numchunks_doppler].reshape((-1, chunksize))
    xchunks_doppler = time_axis_doppler[:chunksize*numchunks_doppler].reshape((-1, chunksize))
    ycenters = ychunks_doppler.mean(axis=1)
    xcenters_doppler = xchunks_doppler.mean(axis=1)
    doppler_average+= ycenters
    
    numchunks_heat = heat_energy_x.size // chunksize 
    ychunks_heat = heat_energy_x[:chunksize*numchunks_heat].reshape((-1, chunksize))
    xchunks_heat = time_axis_heat[:chunksize*numchunks_heat].reshape((-1, chunksize))
    ycenters_heat = xchunks_heat.mean(axis=1)
    xcenters_heat = ychunks_heat.mean(axis=1)
    heating_average+= xcenters_heat       
    
      
    # Calculate the max, min, and means of chunksize-element chunks...
#     max_env = ychunks.max(axis=1)
#     min_env = ychunks.min(axis=1)

    
    # Now plot the bounds and the mean...
#     pyplot.fill_between(xcenters, min_env, max_env, color='gray', 
#                      edgecolor='none', alpha=0.5)
#     pyplot.plot(xcenters, ycenters)
     
#     print np.average(energy_x)

pyplot.semilogy(xcenters_doppler, doppler_average / total_to_average, 'blue', label = 'doppler cooling')
pyplot.semilogy(ycenters_heat, heating_average / total_to_average, 'red', label = 'doppler heating')
pyplot.grid(True, which = 'both')
pyplot.legend()
pyplot.show()