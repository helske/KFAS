\name{approxSSM}
\alias{approxSSM}
\title{Linear Gaussian Approximation for Exponential Family State Space Model}
\usage{
approxSSM(model, theta, maxiter = 50, tol = 1e-15)
}
\arguments{
  \item{model}{A non-Gaussian state space model object of
  class \code{SSModel}.}

  \item{theta}{Initial values for conditional mode theta.}

  \item{maxiter}{The maximum number of iterations used in
  approximation Default is 50.}

  \item{tol}{Tolerance parameter for convergence checks.
  Iterations are continued until
  \eqn{tol>abs(dev_{old}-dev_{new})/(abs(dev_{new})+0.1))}.}
}
\value{
An object which contains the approximating Gaussian state
space model with following additional components:
\item{thetahat}{Mode of \eqn{p(\theta|y)}. }
\item{iterations}{Number of iterations used. }
}
\description{
Function \code{approxSMM} computes the linear Gaussian
approximation of a state space model where observations
follow an exponential family distribution.
}
\details{
The linear Gaussian approximating model is defined by
\deqn{\tilde y_t = Z_t \alpha_t + \epsilon_t, \quad
\epsilon_t \sim N(0,\tilde H_t),}{ytilde[t] = Z[t]\alpha[t]
+ \epsilon[t], \epsilon[t] ~ N(0,Htilde[t]),}
\deqn{\alpha_{t+1} = T_t \alpha_t + R_t \eta_t, \quad
\eta_t \sim N(0,Q_t),}{\alpha[t+1] = T[t]\alpha[t] +
R[t]\eta[t], \eta[t] ~ N(0,Q[t]),} and \eqn{\alpha_1 \sim
N(a_1,P_1)}{\alpha[1] ~ N(a[1],P[1])}, where \eqn{\tilde
y}{ytilde} and \eqn{\tilde H}{Htilde} are chosen in a way
that the linear Gaussian approximating model has the same
conditional mode of \eqn{\theta=Z\alpha} given the
observations \eqn{y} as the original non-gaussian model.
Models also have a same curvature at the mode.

The approximation of the exponential family state space
model is based on iterative weighted least squares method,
see McCullagh and Nelder (1983) p.31 and Durbin Koopman
(2012) p. 243.
}
\seealso{
Importance sampling of non-Gaussian state space models
\code{\link{importanceSSM}}, construct a \code{SSModel}
object \code{\link{SSModel}}, and examples in
\code{\link{KFAS}}.
}

