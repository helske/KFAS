[![R-CMD-check](https://github.com/helske/KFAS/workflows/R-CMD-check/badge.svg)](https://github.com/helske/KFAS/actions)
[![cran version](http://www.r-pkg.org/badges/version/KFAS)](http://cran.r-project.org/package=KFAS)
[![downloads](http://cranlogs.r-pkg.org/badges/KFAS)](http://cranlogs.r-pkg.org/badges/KFAS)

KFAS: R Package for Exponential Family State Space Models
==========================================================================

Package KFAS provides tools for modelling exponential family state space models such as
structural time series, ARIMA models, generalized linear models and generalized linear mixed models.

[Paper at JSS](https://www.jstatsoft.org/article/view/v078i10)

If you use KFAS in your paper, please cite it properly, see `citation("KFAS")` in R, or above link to the paper.

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

KFAS is available at [CRAN](http://cran.r-project.org/web/packages/KFAS/index.html). You can install the latest development version from the GitHub using the devtools package:

```R
install.packages("devtools")
library(devtools)
install_github("helske/KFAS")
```

See
* help(KFAS) in R for examples, and many more examples under different functions, as well as the [Paper at JSS](https://www.jstatsoft.org/article/view/v078i10).
* [ChangeLog](https://github.com/helske/KFAS/blob/master/ChangeLog) for changes.

