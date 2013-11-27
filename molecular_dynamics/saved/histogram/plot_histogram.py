import numpy as np 
from matplotlib import pyplot


# energies = np.load('energy_left.npy')
energies = np.load('energy_left2.npy')

# pyplot.figure()
# pyplot.plot(energies)


# binned,energies,stuff = pyplot.hist(energies, 100)
# print binned.size
# print energies.size

# pyplot.figure()
# pyplot.plot(energies[1:], binned)
pyplot.hist(energies, 100)

pyplot.show()