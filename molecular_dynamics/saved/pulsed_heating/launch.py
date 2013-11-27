'''
python setup.py build_ext --inplace
'''
from simulation import simulator
import numpy as np
from simulation_parameters import simulation_parameters
from equilbrium_positions import equilibrium_positions as equil
from matplotlib import pyplot

doppler_average = np.zeros(599)
pulsed_heating_averged = np.zeros(599)
doppler_heating_averaged = np.zeros(599)
doppler_heating_averagaed_alt = np.zeros(599)


total_to_average = 5
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
    #now do laser heating with pulsing heating
    p.laser_detuning = 0 * p.transition_gamma
    p.saturation = 2.0
    p.simulation_duration = 200e-6
    p.laser_direction = np.array([1.0, 0.0, 0.0])
    p.pulsed_laser = True
    starting_positions = doppler_positions[:,-1,:]
    starting_velocities =  doppler_positions[:,-1,:] -  doppler_positions[:,-2,:]
    heating_simulator = simulator(p)
    heating_positions, heating_excitations = heating_simulator.simulation(starting_positions, starting_velocities, random_seeding = i)
    #now we do laser heating with blue heating
    p.laser_detuning = 0.5 * p.transition_gamma
    p.saturation = 1.0
    p.simulation_duration = 200e-6
    p.laser_direction = np.array([1.0, 0.0, 0.0])
    p.pulsed_laser = False
    starting_positions = doppler_positions[:,-1,:]
    starting_velocities =  doppler_positions[:,-1,:] -  doppler_positions[:,-2,:]
    doppler_heating_simulator = simulator(p)
    doppler_heating_positions, doppler_heating_excitations = doppler_heating_simulator.simulation(starting_positions, starting_velocities, random_seeding = i)
    #now we do laser heating with blue detuning after pulsing for some time
    p.laser_detuning = 0.5 * p.transition_gamma
    p.saturation = 1.0
    p.simulation_duration = 200e-6
    p.laser_direction = np.array([1.0, 0.0, 0.0])
    starting_positions = heating_positions[:,heating_positions.shape[1] // 4,:]
    starting_velocities =  heating_positions[:,heating_positions.shape[1]//4,:] - heating_positions[:,heating_positions.shape[1]//4-1,:]
    doppler_heating_simulator_alt = simulator(p)
    doppler_heating_positions_alt, doppler_heating_excitations_alt = doppler_heating_simulator_alt.simulation(starting_positions, starting_velocities, random_seeding = i)
    
    
    doppler_velocities = np.diff(doppler_positions, axis = 1) / p.timestep
    heat_velocities = np.diff(heating_positions, axis = 1) / p.timestep
    doppler_heat_velocities = np.diff(doppler_heating_positions, axis = 1) / p.timestep
    doppler_heat_velocoties_alt = np.diff(doppler_heating_positions_alt, axis=1) / p.timestep
    
    time_axis_doppler = np.arange(doppler_positions.shape[1]) * p.timestep * 10**6
    time_axis_heat = time_axis_doppler[-1] + np.arange(heating_positions.shape[1]) * p.timestep * 10**6
    time_axis_heat_alt =  time_axis_doppler[-1] + np.arange(heating_positions.shape[1]) * p.timestep * 10**6
    doppler_energy_x = doppler_velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
    heat_energy_x = heat_velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
    doppler_heat_energy_x = doppler_heat_velocities[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
    doppler_heat_alt_energy_x = doppler_heat_velocoties_alt[0, :, 0]**2 * p.mass * .5 / p.hbar / (2 * np.pi * p.f_x)
    
    pyplot.xlabel(r'$\mu s$', fontsize = 30)
    pyplot.ylabel(r'Motional Quanta x', fontsize = 30)
    pyplot.tick_params(axis='both', labelsize = 20)

    def downsample(x, y, points):
        num = x.size // points 
        x = x[:points*num].reshape((-1, points))
        x = x.mean(axis=1)
        y = y[:points*num].reshape((-1, points))
        y = y.mean(axis=1)
        return x,y
    
    time_doppler_cooling, averaged_energy = downsample(time_axis_doppler[1:], doppler_energy_x, 1000)
    doppler_average +=averaged_energy
    time_pulsed_heating,averaged_energy = downsample(time_axis_heat[1:], heat_energy_x, 1000)
    pulsed_heating_averged+=averaged_energy
    time_doppler_heating_immediate,averaged_energy = downsample(time_axis_heat[1:], doppler_heat_energy_x, 1000)
    doppler_heating_averaged+=averaged_energy
    time_doppler_heating_alt,averaged_energy = downsample(time_axis_heat[1:], doppler_heat_alt_energy_x, 1000)
    time_doppler_heating_alt += 50
    doppler_heating_averagaed_alt+=averaged_energy
    

pyplot.plot(time_doppler_cooling, doppler_average / total_to_average, 'b', label = 'doppler cooling')
pyplot.plot(time_pulsed_heating, pulsed_heating_averged / total_to_average, 'r', label = 'pulsed heating')
pyplot.plot(time_pulsed_heating, doppler_heating_averaged / total_to_average, 'r--', label = 'doppler heating')
pyplot.plot(time_doppler_heating_alt, doppler_heating_averagaed_alt / total_to_average, 'r--', label = 'doppler heating alt')

pyplot.grid(True, which = 'both')
pyplot.legend(loc)
pyplot.show()