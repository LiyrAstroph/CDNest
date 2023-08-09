.. CDNest documentation master file, created by
   sphinx-quickstart on Wed Feb  6 21:19:38 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CDNest's documentation!
==================================

**CDNest** is a C version of diffusive nested sampling proposed by Brendon Brewer (https://github.com/eggplantbren/DNest3) with
minor modifications.

 * It provides a C library **libdnest.so**. Examples for using this library are 
   shown with **model1.c**, **model2.c**, and **model3.c** in the subfolder **src/**.
 * It also provides a python callable module **cydnest** wrapped using Cython. 
   Examples for using **cydnest** in Python are shown with **gauss.py**, **mulgauss.py**, 
   and other Python scripts in the subfolder **tests/**.

.. image:: https://www.zenodo.org/badge/DOI/10.5281/zenodo.2527924.svg
   :target: https://doi.org/10.5281/zenodo.2527924

Documentation
-------------

.. toctree::
   :maxdepth: 1
   
   getting_started.rst
   usage.rst
   tests.rst


Attribution
-----------

If you make use of this code, please cite it as:

.. code-block:: tex

  @misc{yan_rong_li_2018_2527924,
    author       = {Yan-Rong Li},
    title        = {CDNest: A diffusive nested sampling code in C},
    month        = dec,
    year         = 2018,
    doi          = {10.5281/zenodo.2527924},
    url          = {https://doi.org/10.5281/zenodo.2527924}
  }


Authors & License
-----------------

Copyright 2019 Yan-Rong Li.

Licensed under the `MIT License <https://opensource.org/licenses/MIT>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
