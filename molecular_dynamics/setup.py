from distutils.core import setup
from Cython.Build import cythonize
 
 
setup(
  name = 'Simulation',
  ext_modules = cythonize("*.pyx")
)