[![Build Status](https://travis-ci.org/helske/KFAS.png?branch=master)](https://travis-ci.org/helske/KFAS)
[![codecov.io](http://codecov.io/github/helske/KFAS/coverage.svg?branch=master)](http://codecov.io/github/helske/KFAS?branch=master)
[![downloads](http://cranlogs.r-pkg.org/badges/KFAS)](http://cranlogs.r-pkg.org/badges/KFAS)
[![cran version](http://www.r-pkg.org/badges/version/KFAS)](http://cran.r-project.org/package=KFAS)
[![Research software impact](http://depsy.org/api/package/cran/KFAS/badge.svg)](http://depsy.org/package/r/KFAS)

KFAS: R Package for Exponential Family State Space Models
==========================================================================

Package KFAS provides tools for modelling exponential family state space models such as
structural time series, ARIMA models, generalized linear models and generalized linear mixed models.

[Paper at JSS](https://www.jstatsoft.org/article/view/v078i10)

If you use KFAS in your paper, please cite properly, see `citation("KFAS")` in R, or above link to the paper. For (old) collection of papers mentioning KFAS, see: [Who uses KFAS](https://rawgit.com/helske/helske.github.io/master/kfas_citedby.html)

Main features
--------------------------------------------------------------------------

- Kalman filtering
- Fixed interval smoothing (Kalman smoothing)
- Simulation smoothing of Gaussian models
- Importance sampling based inference of non-Gaussian models
- Exact diffuse initialization
- Sequential processing
- Multivariate models with mixed distributions

Most of the algorithms are based on book "Time Series Analysis by State Space Methods" and related articles by J. Durbin and S.J. Koopman.

Current version of KFAS in [CRAN](http://cran.r-project.org/web/packages/KFAS/index.html) is 1.2.8. You can install the latest development version from github using devtools package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/KFAS")
```

See
* help(KFAS) in R for examples
* [ChangeLog](https://github.com/helske/KFAS/blob/master/ChangeLog) for upcoming, already completed changes

