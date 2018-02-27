(**To read this papge better, please use google chrome with the extension "GitHub with MathJax" **)

# DNest_C

C version of diffusive nested sampling proposed by Brendon Brewer ( https://github.com/eggplantbren/DNest3 ).

It proivdes a library libdnest.so. Examples for using this library are shown in model1.c, model2., and model3.c.

# Test 1
The likelihood function is (Brewer et al. 2009, arXiv:0912.2380):
$$ L(x_1, x_2, ..., x_20) = \prod_i^{20}\frac{1}{\sqrt{2\pi v^2}}\exp\left(-\frac{x_i^2}{2v^2}\right) + 100 \prod_i^{20}\frac{1}{\sqrt{2\pi u^2}}\exp\left(-\frac{(x_i-0.031)^2}{2u^2}\right),$$
where $v=0.1$ and $u=0.01$. The true value of evidence is $\log(101)\approx4.6151$. The obtained value by ``DNest_c`` for my running is about 4.5807. Note that ``DNest_c`` is not easy to estimate the uncertainty. Different runnings may give slightly different values. A possible way for estimating the uncertainty of evidence is running ``DNest_C`` many times.

# Test 2
A linear regression (Figure 4 in Brewer et al. 2016,  arXiv:1606.0375).
The obtained evidence by ``DNest_C`` is -175.45708.

# Test 3
Gaussian shells likelihood (Feroz et al. 2008, arXiv:0809.3437):
$$L(\boldsymbol{\theta})= \frac{1}{\sqrt{2\pi w_1^2}}\exp\left[-\frac{(|\boldsymbol{\theta-c_1}|-r_1)^2}{2w_1^2}\right]+\frac{1}{\sqrt{2\pi w_2^2}}\exp\left[-\frac{(|\boldsymbol{\theta-c_2}|-r_2)^2}{2w_2^2}\right],$$
where $w_1=w_2=0.1$, $r_1=r_2=2$, and $\boldsymbol{c_1}=(3, 0)$ and $\boldsymbol{c_2}=(-3, 0)$.
The true value of evidence is -1.75. The obtained value by ``DNest_C`` is about -1.7479.

![Gaussian shells likelihood](https://github.com/liyropt/MyGithubPic/blob/master/dnest_test3.jpg)

# Author
Yan-Rong Li
liyanrong@mail.ihep.ac.cn
Please contact me if there is any problem regarding using this library.