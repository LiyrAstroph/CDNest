import os
from setuptools import setup
from setuptools.extension import Extension
from distutils.command.build import build
from setuptools.command.build_ext import build_ext
import subprocess
import numpy

os.environ["CC"] = "mpicc"

include_dirs = ["./", numpy.get_include(),]
library_dirs = ["./", ]

if os.name == 'nt':  # Windows, assumming MSVC compiler
  libraries = ['dnest']
  compiler_args = ['/Ox', '/fp:fast']
elif os.name == 'posix':  # UNIX, assumming GCC compiler
  libraries = ['m', 'c', 'dnest', 'gsl', 'gslcblas']
  compiler_args = ['-O3', '-ffast-math']


try:
  from Cython.Build import cythonize
except ImportError:
  raise RuntimeError('Cython not found.')

class Build(build):
  def run(self):
    command = "cd ./"
    command += " && make"
    process = subprocess.Popen(command, shell=True)
    process.wait()
    build.run(self)

class BuildExt(build_ext):
  def run(self):
    command = "cd ./"
    command += " && make"
    process = subprocess.Popen(command, shell=True)
    process.wait()
    build_ext.run(self)

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
  cmdclass={'build_ext': BuildExt, 'build':Build},
  setup_requires=['numpy'],
  install_requires=['numpy']
	)