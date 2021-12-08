import os
from setuptools import setup
from setuptools.extension import Extension
from distutils.command.build import build
from distutils.command.clean import clean
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

import subprocess
import numpy
import shutil
from glob import glob
import pkgconfig

# if CC is not set, use the default value
if not os.environ.get("CC"):
  os.environ["CC"] = "mpicc"

def configure_mpi():
  """
  get configurations of mpi
  """
  if pkgconfig.exists('mpich'):
    mpiconf = pkgconfig.parse('mpich')
  elif pkgconfig.exists('ompi'):
    mpiconf = pkgconfig.parse('ompi')
  else:
    raise SystemError("Not found MPICH or OpenMPI installed.")

  return mpiconf

mpiconf = configure_mpi()

basedir = os.path.dirname(os.path.abspath(__file__))
homedir = os.environ['HOME']
include_dirs = [basedir, os.path.join(basedir, "src"), numpy.get_include(),] + mpiconf['include_dirs']
library_dirs = [basedir] + mpiconf['library_dirs']

if os.name == 'nt':  # Windows, assumming MSVC compiler
  libraries = ['dnest']
  compiler_args = ['/Ox', '/fp:fast']
  link_args = []
elif os.name == 'posix':  # UNIX, assumming GCC compiler
  libraries = ['m', 'c', 'gsl', 'gslcblas',] + mpiconf['libraries']
  compiler_args = ['-O3', '-ffast-math'] 
  link_args = []

try:
  from Cython.Build import cythonize
except ImportError:
  raise RuntimeError('Cython not found.')

class Build(build):
  def run(self):
    #command = "cd ./"
    #command += " && make"
    #process = subprocess.Popen(command, shell=True)
    #process.wait()
    build.run(self)

class Install(install):
  def run(self):
    #command = "cd ./"
    #command += " && make"
    #process = subprocess.Popen(command, shell=True)
    #process.wait()
    install.run(self)

class BuildExt(build_ext):
  def run(self):
    #command = "cd ./"
    #command += " && make"
    #process = subprocess.Popen(command, shell=True)
    #process.wait()
    build_ext.run(self)

class Clean(clean):
  '''
  Subclass to remove any files created in an inplace build.
  '''
  def run(self):
    clean.run(self)
    # Clean any build or dist directory
    if os.path.isdir("build"):
      shutil.rmtree("build", ignore_errors=True)
    if os.path.isdir("dist"):
      shutil.rmtree("dist", ignore_errors=True)

src = [os.path.join(basedir, "python", "cydnest", "cydnest.pyx")] + glob(os.path.join(basedir, "src", "dnest*.c"))
headerfiles = [os.path.join(basedir, "python", "cydnest", "cydnest.pyd")] + glob(os.path.join(basedir, "src", "dnest*.h"))
extensions = cythonize([
  Extension("cydnest.cydnest", 
	  sources=src,
    depends=headerfiles,
	  extra_compile_args=compiler_args,
    extra_link_args=link_args,
    include_dirs=include_dirs,
    libraries=libraries,
    library_dirs=library_dirs
    ),
  ], annotate=False)

setup(
	name="cydnest",
  version="0.2.0",
	packages=["cydnest"],
  package_dir={"":"python"},
	ext_modules = extensions,
  description = 'C version of Diffusive Nested Sampling (DNest4) by Brendon J. Brewer',
  author = 'Yan-Rong Li',
  author_email = 'liyanrong@mail.ihep.ac.cn',
  cmdclass={'build_ext': BuildExt, 'build':Build, 'install':Install, 'clean':Clean},
  setup_requires=['numpy', 'mpi4py', 'pkgconfig'],
  install_requires=['numpy', 'mpi4py', 'pkgconfig'],
  license="GSL",
  # install header and library to ~/.local/lib
  #data_files=[(os.path.join(homedir, ".local/lib/"), [basedir+"/libdnest.so"]),
  #            (os.path.join(homedir, ".local/include/"), [basedir+"/dnest.h"]),]
	)