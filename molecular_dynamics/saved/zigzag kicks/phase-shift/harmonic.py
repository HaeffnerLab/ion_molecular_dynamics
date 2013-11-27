import numpy as np
from matplotlib import pyplot

w1 = 1
w2 = (1+3.2e-2)*w1
t = np.linspace(0, 1000, 10000)
pyplot.plot(t, np.cos(w1*t) + np.cos(w2*t))
pyplot.plot(t, np.cos(w1*t) - np.cos(w2*t))
pyplot.show()