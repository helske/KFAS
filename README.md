[![Build Status](https://travis-ci.org/helske/KFAS.png?branch=master)](https://travis-ci.org/helske/KFAS)
[![Coverage Status](https://coveralls.io/repos/helske/KFAS/badge.svg?branch=master)](https://coveralls.io/r/helske/KFAS?branch=master)
[![codecov.io](https://codecov.io/github/helske/KFAS/coverage.svg?branch=master)](https://codecov.io/github/helske/KFAS?branch=master)
[![downloads](http://cranlogs.r-pkg.org/badges/KFAS)](http://cranlogs.r-pkg.org/badges/KFAS)
[![cran version](http://www.r-pkg.org/badges/version/KFAS)](http://cran.r-project.org/package=KFAS)

KFAS: R Package for Exponential Family State Space Models
==========================================================================

Package KFAS provides tools for modelling exponential family state space models such as
structural time series, ARIMA models, generalized linear models and generalized linear mixed models.

[Vignette at CRAN](http://cran.r-project.org/web/packages/KFAS/vignettes/KFAS.pdf)

Main features
--------------------------------------------------------------------------

- Kalman filtering
- Fixed interval smoothing (Kalman smoothing)
- Simulation smoothing of Gaussian models
- Importance sampling of non-Gaussian models
- Exact diffuse initialization
- Sequential processing
- Multivariate models with mixed distributions

Most of the algorithms are based on book "Time Series Analysis by State Space Methods" and related articles by J. Durbin and S.J. Koopman.

Current version of KFAS in [CRAN](http://cran.r-project.org/web/packages/KFAS/index.html) is 1.1.1. You can install the latest development version from github using devtools package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/KFAS")
```

See
* help(KFAS) in R for examples
* [ChangeLog](https://github.com/helske/KFAS/blob/master/ChangeLog) for upcoming, already completed changes
* [TODO](https://github.com/helske/KFAS/blob/master/TODO) for future changes
