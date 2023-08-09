************
Installation
************

The dependent third-party packages includes

  * MPICH: a MPI implementation library

  * GSL: the GNU Scientific library

For python interface, additional packages required:
  
  * Cython

  * MPI4PY

  * Numpy

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

Several Python scripts are provided in the tests subdirectory to illustrate its usage. 