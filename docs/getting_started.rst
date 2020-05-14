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

To create the C library **libnest.so**, edit the paths for library and header files in Makefile and then compile using the following terminal command

.. code-block:: bash

  make

Then add the path of **libnest.so** to the system environment variable "LD_LIBRARY_PATH". Edit "bashrc" file in home directory 
and add a line as 

.. code-block:: bash

  export LD_LIBRARY_PATH=/path/to/libnest.so/:$LD_LIBRARY_PATH


After creating ``libnest.so``, to create the Python module **cydnest.so**, use the terminal command

.. code-block:: python 
  
  python setup.py build_ext --inplace

and add the path to **libnest.so** to your python environment settings,

.. code-block:: bash

  export PYTHONPATH=/path/to/CDNest:$PYTHONPATH

If one want to install **cydnest** to the standard python path, use the command 

.. code-block:: python 
  
  python setup.py install

This command by default assume that **mpicc** is located in the standard path, otherwise, use the command 

.. code-block:: python 
  
  CC=/path/to/mpicc/ python setup.py install

.. note::
  Both **python2** and **python3** are supported.
