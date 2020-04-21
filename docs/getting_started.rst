************
Installation
************

The **MPI libaray** (e.g., MPICH) and **GNU Scientific Library** (GSL; http://www.gnu.org/software/gsl) are required. One needs to install these libraries in advance.
If one wants to use python wrapper, python module **mpi4py** is also required.

To create the C library ``libnest.so``, edit the paths for library and header files in Makefile and then compile using the following terminal command

.. code-block:: bash

  make

Then add the path of ``libnest.so`` to the system environment variable "LD_LIBRARY_PATH". Edit "bashrc" file in home directory 
and add a line as 

.. code-block:: bash

  export LD_LIBRARY_PATH=/path/to/libnest.so/:$LD_LIBRARY_PATH


After creating ``libnest.so``, to create the Python module ``cydnest.so``, use terminal command

.. code-block:: python 
  
  python setup.py buid_ext --inplace

or 

.. code-block:: python 

  python setup.py install

.. note::
  Both **python2** and **python3** are supported.
