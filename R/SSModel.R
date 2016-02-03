#' Create a State Space Model Object of Class SSModel
#'
#' Function \code{SSModel} creates a state space object object of class
#' \code{SSModel} which can be used as an input object for various functions of
#' \code{KFAS} package.
#'
#' Formula of the model can contain the usual regression part and additional
#' functions defining different types of components of the model, named as
#' \code{SSMarima}, \code{SSMcustom}, \code{SSMcycle}, \code{SSMregression},
#' \code{SSMseasonal} and \code{SSMtrend}.
#'
#' For more details, see package vignette (the mathematical notation is somewhat non-readable in ASCII).
#'
#' @export
#' @importFrom stats terms update.formula drop.terms model.response model.matrix delete.response update.formula
#' @rdname SSModel
#' @name SSModel
#' @seealso \code{\link{KFAS}} for examples.
#' @param formula An object of class \code{\link{formula}} containing the
#'   symbolic description of the model. The intercept term can be removed with
#'   \code{-1} as in \code{lm}. In case of trend or differenced arima component the
#'   intercept is removed automatically in order to keep the model identifiable.
#'   See package vignette and examples in \code{\link{KFAS}} for special functions
#'   used in model construction.
#' @param data An optional data frame, list or environment containing the
#'   variables in the model.
#' @param H Covariance matrix or array of disturbance terms
#'   \eqn{\epsilon_t}{\epsilon[t]} of observation equation. Defaults to an identity matrix.
#'   Omitted in case of non-Gaussian distributions (augment the state vector if you want to add
#'   additional noise).
#' @param u Additional parameters for non-Gaussian models. See details in
#'   \code{\link{KFAS}}.
#' @param distribution A vector of distributions of the observations. Default is
#'   \code{rep("gaussian", p)}, where \code{p} is the number of series.
#' @param tol A tolerance parameter used in checking whether \code{Finf} or \code{F} is numerically zero.
#'   Defaults to \code{.Machine$double.eps^0.5}. If smoothing gives negative variances for
#'   smoothed states, try adjusting this.
#' @param index A vector indicating for which series the corresponding
#'   components are constructed.
#' @param type For cycle, seasonal, trend and regression components, character
#'   string defining if \code{"distinct"} or \code{"common"} states are used for
#'   different series.
#' @param Q For arima, cycle and seasonal component, a \eqn{p \times p}{p x p}
#'   covariance matrix of the disturbances (or in the time varying case \eqn{p
#'   \times p \times n}{p x p x n} array), where where $p$ = \code{length(index)}.
#'   For trend component, list of length \code{degree} containing the \eqn{p
#'   \times p} or \eqn{p \times p \times n} covariance matrices. For a custom
#'   component, arbitrary covariance matrix or array of disturbance terms
#'   \eqn{\eta_t}{\eta[t]}
#' @param a1 Optional \eqn{m \times 1}{m x 1} matrix giving the expected value
#'   of the initial state vector \eqn{\alpha_1}{\alpha[1]}.
#' @param P1 Optional \eqn{m \times m}{m x m} matrix giving the covariance
#'   matrix of \eqn{\alpha_1}{\alpha[1]}.  In the diffuse case the non-diffuse
#'   part of \eqn{P_1}{P[1]}.
#' @param P1inf Optional \eqn{m \times m}{m x m} matrix giving the diffuse part
#'   of \eqn{P_1}{P[1]}. Diagonal matrix with ones on diagonal elements which
#'   correspond to the unknown initial states.
#' @param R For a custom and regression components, optional \eqn{m \times k}
#'   system matrix or array of transition equation.
#' @param ar For arima component, a numeric vector containing the autoregressive
#'   coeffients.
#' @param ma For arima component, a numericvector containing the moving average
#'   coeffients.
#' @param d For arima component, a degree of differencing.
#' @param stationary For arima component, logical value indicating whether a
#'   stationarity of the arima part is assumed. Defaults to TRUE.
#' @param Z For a custom component, system matrix or array of observation
#'   equation.
#' @param T For a custom component, system matrix or array of transition
#'   equation.
#' @param period For a cycle and seasonal components, the length of the
#'   cycle/seasonal pattern.
#' @param sea.type For seasonal component, character string defining whether to
#'   use \code{"dummy"} or \code{"trigonometric"} form of the seasonal
#'   component.
#' @param degree For trend component, integer defining the degree of the
#'   polynomial trend. 1 corresponds to local level, 2 for local linear trend
#'   and so forth.
#' @param rformula For regression component, right hand side formula or list of
#'   of such formulas defining the custom regression part.
#' @param remove.intercept Remove intercept term from regression model. Default
#'   is \code{TRUE}. This tries to ensure that there are no extra intercept
#'   terms in the model.
#' @param n Length of the series, only used internally for dimensionality check.
#' @param ynames names of the times series, used internally.
#'
#' @return Object of class \code{SSModel}, which is a list with the following
#'   components:
#'   \item{y}{A n x p matrix containing the observations. }
#'   \item{Z}{A p x m x 1 or p x m x n array corresponding to the system matrix
#'   of observation equation. }
#'   \item{H}{A p x p x 1 or p x p x n array
#'   corresponding to the covariance matrix of observational disturbances
#'   epsilon. }
#'   \item{T}{A m x m x 1 or m x m x n array corresponding to the
#'   first system matrix of state equation. }
#'   \item{R}{A m x k x 1 or m x k x n
#'   array corresponding to the second system matrix of state equation. }
#'   \item{Q}{A k x k x 1 or k x k x n array corresponding to the covariance
#'   matrix of state disturbances eta }
#'   \item{a1}{A m x 1 matrix containing the
#'   expected values of the initial states. }
#'   \item{P1}{A m x m matrix
#'   containing the covariance matrix of the nondiffuse part of the initial
#'   state vector. }
#'   \item{P1inf}{A m x m matrix containing the covariance
#'   matrix of the diffuse part of the initial state vector. 
#'   If \code{P1[i,i]} is non-zero then \code{P1inf[i,i]} is automatically set to zero. }
#'   \item{u}{A n x p
#'   matrix of an additional parameters in case of non-Gaussian model.}
#'   \item{distribution}{A vector of length p giving the distributions of the
#'   observations. }
#'   \item{tol}{A tolerance parameter for the diffuse phase. }
#'   \item{call}{Original call to the function. } In addition, object of class
#'   \code{SSModel} contains following attributes:
#'   \item{names}{Names of the
#'   list components. }
#'   \item{p, m, k, n}{Integer valued scalars defining the
#'   dimensions of the model components. }
#'   \item{state_types}{Types of the
#'   states in the model. }
#'   \item{eta_types}{Types of the
#'   state disturbances in the model. }
#'   \item{tv}{Integer vector stating whether \code{Z},\code{H},\code{T},\code{R} or \code{Q} is
#'    time-varying (indicated by 1 in \code{tv} and 0 otherwise).
#'    If you manually change the dimensions of the matrices you must change this attribute also.}
#' @seealso \code{\link{KFAS}} for examples.
#' @examples
#'
#' # example of using data argument
#' y <- x <- rep(1, 3)
#' data1 <- data.frame(x = rep(2, 3))
#' data2 <- data.frame(x = rep(3, 3))
#'
#' f <- formula(~ -1 + x)
#' # With data missing the environment of formula is checked,
#' # and if not found in there a calling environment via parent.frame is checked.
#'
#' c(SSModel(y ~ -1 + x)["Z"]) # 1
#' c(SSModel(y ~ -1 + x, data = data1)["Z"]) # 2
#'
#' c(SSModel(y ~ -1 + SSMregression(~ -1 + x))["Z"]) # 1
#' c(SSModel(y ~ -1 + SSMregression(~ -1 + x, data = data1))["Z"]) # 2
#' c(SSModel(y ~ -1 + SSMregression(~ -1 + x), data = data1)["Z"]) # 2
#' SSModel(y ~ -1 + x + SSMregression(~ -1 + x, data = data1))["Z"] # 1 and 2
#' SSModel(y ~ -1 + x + SSMregression(~ -1 + x), data = data1)["Z"] # both are 2
#' SSModel(y ~ -1 + x + SSMregression(~ -1 + x, data = data1), data = data2)["Z"] # 3 and 2
#'
#' SSModel(y ~ -1 + x + SSMregression(f))["Z"] # 1 and 1
#' SSModel(y ~ -1 + x + SSMregression(f), data = data1)["Z"] # 2 and 1
#' SSModel(y ~ -1 + x + SSMregression(f,data = data1))["Z"] # 1 and 2
#'
#' rm(x)
#' c(SSModel(y ~ -1 + SSMregression(f, data = data1))$Z) # 2
#' \dontrun{
#' # This fails as there is no x in the environment of f
#' try(c(SSModel(y ~ -1 + SSMregression(f), data = data1)$Z))
#' }
SSModel <- function(formula, data, H, u, distribution,
  tol = .Machine$double.eps ^ 0.5) {

  if (missing(data)) {
    data <- environment(formula)
    tsp_data <- NULL
  } else {
    tsp_data <- tsp(data)
    data <- as.data.frame(data)
  }

  # Modifying formula object, catching special functions
  mf <- mc <- match.call(expand.dots = FALSE)
  mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  mf[[1L]] <- as.name("model.frame")
  mf$na.action <- as.name("na.pass")
  components <- c("SSMregression", "SSMtrend", "SSMseasonal", "SSMcycle",
    "SSMarima", "SSMcustom")

  all_terms <- terms_out <- terms(formula, specials = components, data = data)
  specials <- attr(all_terms, "specials")
  components <- components[!sapply(specials, is.null)]
  if (length(unlist(specials)) > 0) {
    if (length(attr(all_terms, "term.labels")) == length(unlist(specials))){
      all_terms <- terms(update.formula(all_terms, . ~ . + .emptyx.),
        specials = components)
    }
    drops <- which(attr(all_terms, "term.labels") %in%
        rownames(attr(all_terms, "factors"))[unlist(specials)])
    mf$formula <- formula(drop.terms(all_terms, drops, keep.response = TRUE))
    mf$formula <- update.formula(mf$formula, . ~ . - .emptyx., simplify = TRUE)
  }
  # Remove extra intercept
  if ("SSMtrend" %in% components ||
      ("SSMarima" %in% components &&
          isTRUE(eval(attr(all_terms, "variables")[[specials$SSMarima + 1]]$d,
            envir = data, enclos=parent.frame()) > 0))){
    mf$formula <- update.formula(mf$formula, . ~ . - 1)
    remove.intercept<-TRUE
  } else remove.intercept<-FALSE

  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "numeric")
  mt <- attr(mf, "terms")

  vars <- attr(all_terms, "variables")
  reg_in_formula <- as.integer(dim(model.matrix(mt, mf))[2] > 0)
  specials <- unlist(specials)
  lspecials <- length(specials)
  n_blocks <- lspecials + reg_in_formula
  blocks <- vector("list", n_blocks)

  # building y
  if (is.array(y)) {
    p <- dim(y)[2]
    n <- dim(y)[1]
  } else {
    y <- as.array(y)
    if (length(dim(y)) != 2) {
      p <- 1
      n <- length(y)
      dim(y) <- c(n, p)
    }
  }
  y_names <- colnames(y)
  if (is.null(y_names)) {
    y_names <- if (p > 1) {
      colnames(y) <- paste0("y", 1:p)
    } else ""
  }
  class(y) <- if (p > 1)
    c("mts", "ts", "matrix") else "ts"
  if (is.null(tsp(y))) {
    if (!is.null(tsp_data)) {
      tsp(y) <- tsp_data
    } else tsp(y) <- c(1, n, 1)
  }
  # defining the distributions of y
  if (missing(distribution)) {
    distribution <- rep("gaussian", length = p)
  } else {
    if (length(distribution) == 1 | length(distribution) == p) {
      distribution <- rep(pmatch(x = distribution, table = c("gaussian", "poisson",
        "binomial", "gamma", "negative binomial"), duplicates.ok = TRUE),
        length = p)
    } else stop("Length of the argument 'distribution' must be either 1 or p, the number of series.")
    if (any(is.na(distribution)))
      stop("Misspeficied distribution, only 'gaussian', 'poisson', 'binomial', 'gamma', and 'negative binomial' are allowed.")
    distribution <- c("gaussian", "poisson", "binomial", "gamma", "negative binomial")[distribution]
  }
  # building H and u
  if (all(distribution == "gaussian")) {
    if (!missing(H)) {
      if (length(H) == 1)
        dim(H) <- c(1, 1)
      dims <- dim(H)
      if (dims[1] != p || dims[2] != p)
        stop("Misspecified H, argument H must be NULL, a scalar, p x p matrix or p x p x n array, where p is the number of time series.")
      H <- array(H, dim = c(p, p, 1 + (n - 1) * (max(dims[3], 0, na.rm = TRUE) >
          1)))
    } else {
      H <- array(diag(p), dim = c(p, p, 1))
    }
    u <- "Omitted"
  } else {
    if (!missing(H))
      warning("H ignored as model contains non-gaussian series.")
    H <- "Omitted"
    if (missing(u)) {
      u <- array(1, c(n, p))
    } else {
      if (is.data.frame(u)) {
        u <- data.matrix(u)
        if (!identical(dim(u), c(n, p)))
          stop("Mispecified u, argument u must be either vector of length p, or n x p matrix, where p is the number of time series.")
      } else u <- matrix(u, n, p, byrow = is.vector(u))
      storage.mode(u) <- "double"
    }
    class(u) <- class(y)
    tsp(u) <- tsp(y)
  }
  # building model components by calling appropriate component functions
  if (reg_in_formula) {
    blocks[[1]] <- SSMregression(rformula = formula(delete.response(mt)), data = mf,
      index = 1:p, n = n, ynames = if (p > 1)
        y_names,remove.intercept=remove.intercept)
    blocks[[1]]$state_types <- "regression"
  }
  if (lspecials > 0)
    for (i in 1:lspecials) {

      comp <- vars[[1 + specials[i]]]
      comp <- match.call(definition = eval(comp[[1]],envir=data,enclos=parent.frame()), call = comp)
      if (is.null(comp$index)) {
        comp$index <- 1:p
      } else if (!all(eval(comp$index,envir=data,enclos=parent.frame()) %in% (1:p)))
        stop("Index must have values between 1 to p. ")
      comp$n <- n
      if (comp[[1]] != "SSMcustom" && p > 1 &&
          is.null(comp$ynames) && (is.null(comp$type) || comp$type != "common"))
        comp$ynames <- y_names[eval(comp$index,envir=data,enclos=parent.frame())]

      blocks[[i + reg_in_formula]] <- eval(comp, envir=data,enclos=parent.frame())
      # use level and slope instead of generic trend
      if((s_type <- substr(as.character(comp[[1]]),
        start = 4, stop = 15L))=="trend"){
        s_type <- substr(blocks[[i + reg_in_formula]]$state_names, start = 1, stop = 5L)
      }
      blocks[[i + reg_in_formula]]$state_types <- s_type
    }

  # building combined model arrays
  cum_m <- c(0, cumsum(unlist(sapply(blocks, "[", "m"))))
  cum_k <- c(0, cumsum(unlist(sapply(blocks, "[", "k"))))
  m <- max(cum_m)
  k <- max(cum_k, 1)
  tv<-numeric(5)
  tv[1] <- max(0, unlist(sapply(blocks, "[", "tvz")))
  tv[2] <- any(distribution != "gaussian") || (dim(H)[3] > 1)
  tv[3] <- max(0, unlist(sapply(blocks, "[", "tvt")))
  tv[4] <- max(0, unlist(sapply(blocks, "[", "tvr")))
  tv[5] <- max(0, unlist(sapply(blocks, "[", "tvq")))
  Z <- array(0, c(p, m, 1 + (n - 1) * tv[1]))
  T <- array(0, c(m, m, 1 + (n - 1) * tv[3]))
  R <- array(0, c(m, k, 1 + (n - 1) * tv[4]))
  Q <- array(0, c(k, k, 1 + (n - 1) * tv[5]))
  P1 <- P1inf <- matrix(0, m, m)
  a1 <- matrix(0, nrow = m)
  state_names <- unname(unlist(sapply(blocks, "[", "state_names")))
  rownames(a1) <- rownames(T) <- colnames(T) <- colnames(Z) <-
    rownames(R) <- rownames(P1) <- colnames(P1) <- rownames(P1inf) <-
    colnames(P1inf) <- state_names
  rownames(Z) <- colnames(y)
  state_types <- character(m)
  eta_types <- character(k)
  for (i in 1:n_blocks) {
    Z[blocks[[i]]$index, (cum_m[i] + 1):cum_m[i + 1], ] <- blocks[[i]]$Z
    T[(cum_m[i] + 1):cum_m[i + 1], (cum_m[i] + 1):cum_m[i + 1], ] <- blocks[[i]]$T
    R[(cum_m[i] + 1):cum_m[i + 1], seq_len(cum_k[i + 1] - cum_k[i]) + cum_k[i],
      ] <- blocks[[i]]$R
    Q[seq_len(cum_k[i + 1] - cum_k[i]) + cum_k[i], seq_len(cum_k[i + 1] - cum_k[i]) +
        cum_k[i], ] <- blocks[[i]]$Q
    a1[(cum_m[i] + 1):cum_m[i + 1], ] <- blocks[[i]]$a1
    P1[(cum_m[i] + 1):cum_m[i + 1], (cum_m[i] + 1):cum_m[i + 1]] <- blocks[[i]]$P1
    P1inf[(cum_m[i] + 1):cum_m[i + 1], (cum_m[i] + 1):cum_m[i + 1]] <- blocks[[i]]$P1inf
    state_types[(cum_m[i] + 1):cum_m[i + 1]] <- eta_types[seq_len(cum_k[i + 1] -
        cum_k[i]) + cum_k[i]] <- blocks[[i]]$state_types
  }
  if (all(dim(R) == c(1, 1, 1)) && R[1] == 0)
    R[1] <- 1
  model <- list(y = y, Z = Z, H = H, T = T, R = R, Q = Q, a1 = a1, P1 = P1, P1inf = P1inf,
    u = u, distribution = distribution, tol = tol)
  class(model) <- "SSModel"
  attr(model, "p") <- as.integer(p)
  attr(model, "m") <- as.integer(m)
  attr(model, "k") <- as.integer(k)
  attr(model, "n") <- as.integer(n)
  names(tv) <- c("Z","H","T","R","Q")
  attr(model, "tv") <- as.integer(tv)
  attr(model, "state_types") <- state_types
  attr(model, "eta_types") <- eta_types
  model$call <- mc
  model$terms <- terms_out
  invisible(model)
}
