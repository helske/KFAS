[![Build Status](https://travis-ci.org/helske/KFAS.png?branch=master)](https://travis-ci.org/helske/KFAS)
[![Coverage Status](https://coveralls.io/repos/helske/KFAS/badge.svg?branch=master)](https://coveralls.io/r/helske/KFAS?branch=master)

KFAS: R Package for Exponential Family State Space Models
==========================================================================

Package KFAS provides tools for modelling exponential family state space models such as
structural time series, ARIMA models, generalized linear models and generalized linear mixed models.

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

Current version of KFAS in CRAN is 1.0.4-1. You can install the latest development version from github using devtools package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/KFAS")
```

See
* help(KFAS) in R for examples
* [ChangeLog](https://github.com/helske/KFAS/blob/master/ChangeLog) for upcoming, already completed changes
* [TODO](https://github.com/helske/KFAS/blob/master/TODO) for future changes
