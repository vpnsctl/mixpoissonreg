# The *mixpoissonreg* package

### Dev badges:
<!-- badges: start -->
[![](https://img.shields.io/badge/devel%20version-1.0.0-blue.svg)](https://github.com/vpnsctl/mixpoissonreg/)
[![R build status](https://github.com/vpnsctl/mixpoissonreg/workflows/R-CMD-check/badge.svg)](https://github.com/vpnsctl/mixpoissonreg/actions)
[![Codecov test coverage](https://codecov.io/gh/vpnsctl/mixpoissonreg/branch/main/graph/badge.svg?token=JXDPBKWGYE)](https://codecov.io/gh/vpnsctl/mixpoissonreg)
<!-- badges: end -->
---

The *mixpoissonreg* package deals with regression models with response variables being count data. 
It is aimed towards overdispersed data. It currently provides implementation of two regression models: Negative Binomial and Poisson Inverse Gaussian regression models. For both of these models there are two estimation procedures implemented: the EM-algorithm and direct maximization of the 
log-likelihood function.  The definition of both these regression models along with the EM-algorithm approach
for them can be found at [Barreto-Souza and Simas (2016)](https://doi.org/10.1007/s11222-015-9601-6). The direct maximization
of the log-likelihood is only recommended when the EM-algorithm takes too long to converge.
Several global and local influence measures are implemented and easy to interpret through plot visualization. There are plot implementations
using base-R and *ggplot2*. For further details, and also more tutorials, we refer the reader to the **vignettes** whose links can be found below
or by using `vignette(mixpoissonreg)`.

## Installing the *mixpoissonreg* package

The latest **stable** development version can be installed from github:

```{r}
#install.packages("devtools")
devtools::install_github("vpnsctl/mixpoissonreg")
```

## Dependencies

The *mixpoissonreg* package imports the `Rfast` package. The `Rfast` package needs the `GSL` library to be installed. If you are using linux, 
you most likely will need to install the `GSL` library. 

On Ubuntu (version 18.04 or newer), run:
```{bash}
sudo apt install libgsl-dev
```

On RedHat or Fedora, run:
```{bash}
yum install gsl gsl-devel
```

The *mixpoissonreg* package also imports the following packages:
`pbapply`, `Formula`,  `gamlss`, `gamlss.dist`, `rlang`, `statmod`, `lmtest`, `generics`, `magrittr`, `tibble`, `dplyr`, `ggplot2`,
`ggrepel`, `gridExtra`

## Useful Resources

* The [*mixpoissonreg* home page](https://vpnsctl.github.io/mixpoissonreg/) contains the documentation (with examples), vignettes, etc.

* The *mixpoissonreg* Documentation can be accessed directly through [this link](https://vpnsctl.github.io/mixpoissonreg/reference/index.html)

* Currently, the *mixpoissonreg* package has the following vignettes:

    * [Getting started: Introduction to the *mixpoissonreg* package](https://vpnsctl.github.io/mixpoissonreg/articles/mixpoissonreg.html)
    * [Analyzing overdispersed count data with the *mixpoissonreg* package](https://vpnsctl.github.io/mixpoissonreg/articles/tutorial-mixpoissonreg.html)
    * [Maximum-likelihood estimation with the *mixpoissonreg* package](https://vpnsctl.github.io/mixpoissonreg/articles/ml-mixpoissonreg.html)
    * [Building and customizing base-R diagnostic plots with the mixpoissonreg package](https://vpnsctl.github.io/mixpoissonreg/articles/plots-mixpoissonreg.html)
    * [Building and customizing ggplot2-based diagnostic plots with the mixpoissonreg package](https://vpnsctl.github.io/mixpoissonreg/articles/ggplot2-plots-mixpoissonreg.html)
    * [Global and local influence analysis with the *mixpoissonreg* package](https://vpnsctl.github.io/mixpoissonreg/articles/influence-mixpoissonreg.html)
    * [Confidence and prediction intervals with  the *mixpoissonreg* package](https://vpnsctl.github.io/mixpoissonreg/articles/intervals-mixpoissonreg.html)
    * [*mixpoissonreg* in the *tidyverse*](https://vpnsctl.github.io/mixpoissonreg/articles/tidyverse-mixpoissonreg.html)

* Theoretical details can be obtained in [Barreto-Souza and Simas (2016)](https://doi.org/10.1007/s11222-015-9601-6) and also in its [Supplementary Material](https://link.springer.com/article/10.1007%2Fs11222-015-9601-6#Sec23)

## Basic usage

The usage of the *mixpoissonreg* package is analogous to the usage of standard regression functions and packages in *R*:

```{r}
library(mixpoissonreg)
fit <- mixpoissonreg(daysabs ~ prog + math + gender | prog, data = Attendance)
fit
#> 
#> Negative Binomial Regression - Expectation-Maximization Algorithm
#> 
#> Call:
#> mixpoissonreg(formula = daysabs ~ prog + math + gender | prog, 
#>     data = Attendance)
#> 
#> Coefficients modeling the mean (with log link):
#>    (Intercept)   progAcademic progVocational           math     gendermale 
#>    2.752418105   -0.424045918   -1.238680841   -0.006791582   -0.257120287 
#> Coefficients modeling the precision (with log link):
#>    (Intercept)   progAcademic progVocational 
#>       1.109533      -1.083443      -1.517825

summary(fit)
#> 
#> Negative Binomial Regression - Expectation-Maximization Algorithm
#> 
#> Call:  
#> mixpoissonreg(formula = daysabs ~ prog + math + gender | prog, 
#>     data = Attendance)
#> 
#> 
#> Pearson residuals:
#>      RSS      Min       1Q   Median       3Q      Max 
#> 323.6397  -1.0648  -0.7204  -0.3654   0.3064   4.8776 
#> 
#> Coefficients modeling the mean (with  link):
#>                 Estimate Std.error z-value Pr(>|z|)    
#> (Intercept)     2.752418  0.152613  18.035  < 2e-16 ***
#> progAcademic   -0.424046  0.132369  -3.204  0.00136 ** 
#> progVocational -1.238681  0.173155  -7.154 8.45e-13 ***
#> math           -0.006792  0.002274  -2.986  0.00283 ** 
#> gendermale     -0.257120  0.116934  -2.199  0.02789 *  
#> 
#> Coefficients modeling the precision (with  link):
#>                Estimate Std.error z-value Pr(>|z|)    
#> (Intercept)      1.1095    0.2783   3.987 6.68e-05 ***
#> progAcademic    -1.0834    0.3074  -3.524 0.000424 ***
#> progVocational  -1.5178    0.3395  -4.471 7.79e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
#> 
#> Efron's pseudo R-squared:  0.1861217 
#> Number of iterations of the EM algorithm = 1930
```

There are also several methods implemented:
- Visualization: `plot` (base-R visualization) ; `autoplot` (*ggplot2* visualization) ; `local_influence_plot` (base-R visualization) ; `local_influence_autoplot`  (*ggplot2* visualization)
- Inference: `coeftest` ; `coefci` ; `lmtest::lrtest` (works with the default method) ; `lmtest::waldtest` (works with the default method) ; `predict` 
- Residual analysis: `residuals`
- Global influence analysis: `influence` ; `cooks.distance` ; `hatvalues`
- Local influence analysis: `local_influence`
- Tidyverse compatibility: `augment` ; `glance` ; `tidy` ; `tidy_local_influence` ; `local_influence_benchmarks`
