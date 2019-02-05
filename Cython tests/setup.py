try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
# import numpy as np

# from distutils.core import setup
# from distutils.extension import Extension
from Cython.Build import cythonize

setup(ext_modules=cythonize('ising.pyx'))
# include_dirs=[numpy.get_include()]



# setup(
#     ext_modules=[
#         Extension("ising", ["ising.c"],
#                   include_dirs=[np.get_include()]),
#     ],
# )
# setup(
#     ext_modules=cythonize("ising.pyx"),
#     include_dirs=[np.get_include()]
# )
