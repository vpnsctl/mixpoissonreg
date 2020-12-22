# The *mixpoissonreg* package

### Dev badges:
[![](https://img.shields.io/badge/devel%20version-1.0.0-blue.svg)](https://github.com/vpnsctl/mixpoissonreg/main)
[![R build status](https://github.com/vpnsctl/mixpoissonreg/workflows/R-CMD-check/badge.svg)](https://github.com/vpnsctl/mixpoissonreg/actions)
[![Build Status](https://travis-ci.com/vpnsctl/mixpoissonreg.svg?branch=main)](https://travis-ci.com/vpnsctl/mixpoissonreg)
---

The *mixpoissonreg* package deals with regression models with response variables being count data. 
It currently provides implementation of two regression models:  Poisson Inverse Gaussian regression
Negative Binomial. For both of these models, the estimation is primarily
done with the EM-algorithm. The definition of both these regression models along with the EM-algorithm approach
for them can be found in <http://doi.org/10.1007/s11222-015-9601-6>. We also provide maximum-likelihood estimation
for both of these models, which is only recommended when the EM-algorithm takes too long to converge.

## Installing the *mixpoissonreg* package from this repository:

To install the *mixpoissonreg* package from this repository, just run the command:

```{r}
#install.packages("devtools")
devtools::install_github("vpnsctl/mixpoissonreg")
```

To install the *mixpoissonreg* package from this repository **with vignettes**, run the command:
```{r}
#install.packages("devtools")
devtools::install_github("vpnsctl/mixpoissonreg", build_vignettes = TRUE)
```

This repository will always contain the most recent version of the *mixpoissonreg* package.

## Dependencies

The *mixpoissonreg* package depends on the *Rfast* package. The *Rfast* package needs the *GSL* library to be installed. If you are using linux, you most likely will need to install the *GSL* library. 

On Ubuntu (version 18.04 or newer), run:
```{bash}
sudo apt install libgsl-dev
```

On RedHat or Fedora, run:
```{bash}
yum install gsl gsl-devel
```

## Vignettes

The *mixpoissonreg* package currently has the following vignettes:

* [Introduction to the *mixpoissonreg* package](https://rpubs.com/alexandrebsimas/intro-mixpoissonreg)
* [Analyzing overdispersed count data with the *mixpoissonreg* package](https://rpubs.com/alexandrebsimas/tutorial-mixpoissonreg)
* [Maximum-likelihood estimation with the *mixpoissonreg* package](https://rpubs.com/alexandrebsimas/ml-mixpoissonreg)
* [Diagnostic plots with the *mixpoissonreg* package](https://rpubs.com/alexandrebsimas/plots-mixpoissonreg)
* [Global and local influence analysis with the *mixpoissonreg* package](https://rpubs.com/alexandrebsimas/influence-mixpoissonreg)
* [Confidence and predictions intervals with  the *mixpoissonreg* package](https://rpubs.com/alexandrebsimas/intervals-mixpoissonreg)
* [*mixpoissonreg* in the tidyverse](https://rpubs.com/alexandrebsimas/tidyverse-mixpoissonreg)

## Basic usage

The usage of the *mixpoissonreg* package is analogous to the usage of standard regression functions and packages in *R*:

For further details we refer the reader to the **vignettes** whose link can be found above.
