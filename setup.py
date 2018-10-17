import os
from setuptools import setup
from setuptools.extension import Extension

os.environ["CC"] = "mpicc"

include_dirs = ["/home/liyropt/Projects/GIT/DNest",]
library_dirs = ["/home/liyropt/Projects/GIT/DNest",]

if os.name == 'nt':  # Windows, assumming MSVC compiler
  libraries = ['dnest']
  compiler_args = ['/Ox', '/fp:fast']
elif os.name == 'posix':  # UNIX, assumming GCC compiler
  libraries = ['m', 'dnest']
  compiler_args = ['-O3', '-ffast-math']


try:
  from Cython.Build import cythonize
except ImportError:
  raise RuntimeError('Cython not found.')

extensions = cythonize([
  Extension("cydnest", 
	  sources=["cydnest.pyx",],
	  extra_compile_args=compiler_args,
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs
    ),
  ], annotate=False)

setup(
	name="cydnest",
	packages="cydnest",
	ext_modules = extensions,
  description = 'C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer',
  author = 'Yan-Rong Li',
  author_email = 'liyanrong@mail.ihep.ac.cn',
	)