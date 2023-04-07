******
Tests
******

**CDNest** provides several tests to guide its useage in the subdirectory **tests/**.

Test I --- model1.c
===================

The likelihood function is (Brewer et al. 2009, arXiv:0912.2380):

.. math::
  
   L(x_1, x_2, ..., x_20) = \prod_i^{20}\frac{1}{\sqrt{2\pi v^2}}\exp\left(-\frac{x_i^2}{2v^2}\right) + 100 \prod_i^{20}\frac{1}{\sqrt{2\pi u^2}}\exp\left(-\frac{(x_i-0.031)^2}{2u^2}\right),

where :math:`v=0.1` and :math:`u=0.01`. The true value of evidence is :math:`\log(101)\approx4.6151`. The obtained value by CDNest for my running is about 4.5807. Note that CDNest is not easy to estimate the uncertainty. Different runnings may give slightly different values. A possible way for estimating the uncertainty of evidence is running CDNest many times.

To run this test, using the following command:

.. code-block:: bash
  
  mpiexec -n np ./dnest 1

Test II --- model2.c
====================

A linear regression (Figure 4 in Brewer et al. 2016, arXiv:1606.0375). The obtained evidence by CDNest is -175.45708.

To run this test, using the following command:

.. code-block:: bash
  
  mpiexec -n np ./dnest 2

Test III --- model3.c
=====================

Gaussian shells likelihood (Feroz et al. 2008, arXiv:0809.3437):

.. math::
  
  L(\boldsymbol{\theta})= \frac{1}{\sqrt{2\pi w_1^2}}\exp\left[-\frac{(|\boldsymbol{\theta-c_1}|-r_1)^2}{2w_1^2}\right]+\frac{1}{\sqrt{2\pi w_2^2}}\exp\left[-\frac{(|\boldsymbol{\theta-c_2}|-r_2)^2}{2w_2^2}\right],

where :math:`w_1=w_2=0.1, r_1=r_2=2, \boldsymbol{c_1}=(3, 0)` and :math:`\boldsymbol{c_2}=(-3, 0)`. The true value of evidence is -1.75. The obtained value by CDNest is about -1.7479.

.. figure:: _static/fig_test3.jpg
  :scale: 100 %
  :align: center

To run this test, using the following command:

.. code-block:: bash
  
  mpiexec -n np ./dnest 3

Test IV in Python --- gauss.py
================================

Generate samples form a Gaussian distribution.

.. figure:: _static/fig_gau.jpg
  :scale: 100%
  :align: center

To run this test, using the following command:

.. code-block:: bash
  
  mpiexec -n np python gauss.py

Test V in Python --- mulgauss.py
================================

Multidimensional Gaussian.

The true value of evidence is -11.51, and the obtained value by CDNest is -11.52.

To run this test, using the following command:

.. code-block:: bash
  
  mpiexec -n np python mulgauss.py

Test VI in Python --- rastrigin.py
==================================

The two-dimensional Rastrigin test function is defined by 

.. math::
  
  f(\theta) = An + \sum_{i=1}^{n}[\theta_i^2 - A \cos(2\pi\theta_i)],

  A = 10, n=2, \theta_i \sim [-5.12, 5.12]

The likelihood is defined to be :math:`L = \exp[-f(\theta)]`.

.. figure:: _static/fig_rastrigin.jpg
  :scale: 50%
  :align: center

To run this test, using the following command:

.. code-block:: bash
  
  mpiexec -n np python rastrigin.py

Test V in Python --- gauss_plateau.py
=====================================

A 2D clipped Gaussian likelihood, defined as 

.. math::

  &L(x, y) = \frac{1}{2\pi} \exp\left(-\frac{x^2+y^2}{2} \right) ~for~ x^2+y^2 <= 4

  &L(x, y) = \frac{1}{2\pi} \exp\left(- 2 \right) ~for~ x^2+y^2 > 4

We set a uniform prior for both :math:`x` and :math:`y` in a range (-5, 5). The true 
evidence is :math:`\log(Z) = -2 + \log(3/100)\approx-5.51`.  The evidence obtained 
value by CDNest is -5.51.

.. figure:: _static/fig_gauss_plateau.jpg 
  :scale: 30%
  :align: center

To run this test, using the following command:

  .. code-block:: bash
    
    mpiexec -n np python gauss_plateau.py

Test VI in Python --- cauchy.py
===============================

The 48-d Cauchy likelihood is defined as 

.. math::
  &L = \prod_{i=1}^{n} \frac{1}{2}\left[ {\rm Cauchy}(x_i|\mu, \sigma) + {\rm Cauchy}(x_i|-\mu, \sigma)\right] 

  & \mu = 5, \sigma = 1, n=48

We set a uniform prior over a range [-100, 100] for all the parameters. 
The true evidence is :math:`\log(Z) = 48(-\ln(200) + \ln(2\tan^{-1}(100)/\pi)) \approx 254.626`. The obtained evidence 
by CDNest is -254.654. 

.. figure:: _static/fig_cauchy.jpg 
  :scale: 50%
  :align: center

To run this test, using the following command:

  .. code-block:: bash
    
    mpiexec -n np python cauchy.py

Test VII in Python --- diamond_ring.py
======================================

The diamond ring problem, with the likelihood as (see Buchner 2023, arXiv:2101.09675)

.. math:: 
  & L = N_1 + 100 N_2

  & N_i = \frac{1}{\sqrt{2\pi}w_i} \exp \left[-\frac{1}{2}\left(\frac{d_i-r_i}{w_i^2}\right)^2 \right]

  & r_1 = 10^{-11}, w_1 = 0.4\times r_1, d_1 = \sqrt{x^2+y^2}

  & r_2 = r_1/40, w_1 = w_2/40, d_2 = \sqrt{(x+r_2)^2+y^2}

We set a uniform prior over a range [-1, 1] for both :math:`x` and :math:`y`. Note that small scales of the likelihood 
peak regions compared to the priors.

The true evidence is :math:`\log(Z)\approx-23.62`. The obtained evidence by CDNest is -23.89 by one run 
(the value subjecting to fluctuations). 

.. figure:: _static/fig_diamond.jpg 
  :scale: 50%
  :align: center

To run this test, using the following command:

  .. code-block:: bash
    
    mpiexec -n np python diamond_ring.py