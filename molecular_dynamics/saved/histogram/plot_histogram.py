import numpy as np 
from matplotlib import pyplot


energies = np.load('energy_left3.npy')#25 ions, long
# energies = np.load('energy_left5ms_5ions.npy');energies = energies[0] #5 ions long


pyplot.subplot(311)
pyplot.title('5 ions energy dynamics')
pyplot.plot(energies)
pyplot.xlim(0,1000)
pyplot.subplot(312)
pyplot.title('energy histogarm')
pyplot.hist(energies, 100)
pyplot.subplot(313)
pyplot.title('autocorrelation')
pyplot.acorr(energies, maxlags = 1000)
pyplot.xlim(0,1000)
pyplot.tight_layout()
pyplot.show()