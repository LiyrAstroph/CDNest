************
Installation
************

The dependent third-party packages includes

  * MPICH: `an MPI implementation library <https://www.mpich.org/>`_

  * GSL: `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
  
In Linux systems, one can use the package manager (like ``dnf``, which is different among different distributions) 
to install the above packages 

.. code-block:: bash 

  dnf install mpich mpich-devel gsl gsl-devel

In Mac systems, one can use `homebrew <https://brew.sh>`_ package manager to install the above packages 

.. code-block:: bash 

  brew install mpich gsl

**For python interface**, additional packages required:
  
  * Cython

  * MPI4PY

  * Numpy

Using ``pip`` to install these packages  

.. code-block:: bash 

  pip install cython mpi4py numpy


C/C++ library: libdnest.so
==========================
To create the C library **libdnest.so**, edit the paths for library and header files in Makefile and 
then compile using the following terminal command

.. code-block:: bash

  make

Then add the path of **libdnest.so** to the system environment variable "LD_LIBRARY_PATH". 
Edit "bashrc" file in home directory and add a line as 

.. code-block:: bash

  export LD_LIBRARY_PATH=/path/to/libdnest.so/:$LD_LIBRARY_PATH


Python package: cydnest
========================
After creating ``libnest.so``, to create the Python package **cydnest**, use the terminal command

.. code-block:: bash 
  
  python setup.py build_ext --inplace

This will generate Python-callable package **cydnest** in the local path **./build/libXXX/cydnest**, where 
**XXX** depend on the compiling system. (In my Fedora 33 laptop, the path name is lib.linux-x86_64-3.9)
Add this path to your python environment settings,

.. code-block:: bash

  export PYTHONPATH=/path/to/CDNest/build/libXXX/cydnest:$PYTHONPATH

If one wants to install **cydnest** to the standard python path, use the command 

.. code-block:: bash 
  
  sudo python setup.py install

or 

.. code-block:: bash 
  
  python setup.py install --user

The former command generally installs cydnest to the path like **/usr/libXX/pythonXX/site-packages**, 
while the latter to the path like **~/.local/libXX/pythonXX/site-packages/**, where **XX** depends 
on the system environment.

The above commands by default assume that **mpicc** is located in the standard path, otherwise, use the command 

.. code-block:: python 
  
  CC=/path/to/mpicc/ python setup.py install

.. note::
  Both **python2** and **python3** are supported.

To use the package **cydnest**,  import it in a Python scirpt as 

.. code-block:: python 

  import cydnest

Several Python scripts are provided in the **tests** subdirectory to illustrate its usage. 

Outputs
=========

**CDNest** forces writing to the computer's disk every **num_max_saves/50** steps, so one does not need to wait until the end of 
running, but instead, one can inspect the results during the running and adjust the options when necessary.

The default output files of **CDNest** are as follows.

- **sample.txt**
  
  The ancillary file to record saved parameters during the sampling. Each row corresponds to a set of parameters 
  and the number of rows is the **num_max_saves** option. 

- **sample_info.txt**
  
  The level assignment, log likelihood, tiebreaker, id of each set of parameters in **sample.txt**. 
  Here tiebreaker is a random number is used to distinguish parameter sets with the same likelihood,
  and id is a number marking which CUP cores the saved parameters come from. 
  
- **sampler_state.txt**
  
  Record the number of levels and the number of rows of parameters saved.

- **levels.txt**
  
  Record the level informations.

- **limits.txt**
  
  Record the limits of each parameter at each level.

- **posterior_sample.txt**
  
  **The posterior sample of parameters.**

- **posterior_sample_info.txt**
  
  **The likelihod of posterior sample.**
