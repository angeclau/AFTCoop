% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_functions.R
\name{hessian_coop}
\alias{hessian_coop}
\title{Calculate Hessian and Lipschitz Constant for AFT Models}
\usage{
hessian_coop(X, Y, eta, delta, sigma, rho, lambda, coop, n, model)
}
\arguments{
\item{X}{Matrix of covariates}

\item{Y}{Vector of observed times (response variable)}

\item{eta}{Vector of linear predictors}

\item{delta}{Vector of censoring indicators (1 for event, 0 for censored)}

\item{sigma}{Scale parameter}

\item{rho}{Cooperation parameter}

\item{lambda}{Regularization parameter}

\item{coop}{Cooperative matrix}

\item{n}{Sample size}

\item{model}{Character string specifying the distribution ("weibull", "lognormal", or "loglogistic")}
}
\value{
A list containing:
\itemize{
\item L: Lipschitz constant (maximum eigenvalue of Hessian)
}
}
\description{
Computes the Hessian matrix and its maximum eigenvalue (Lipschitz constant)
for the negative log-likelihood plus cooperative penalty.
}
\details{
The function:
\itemize{
\item Computes second derivatives specific to each distribution
\item Includes cooperative penalty term in Hessian
\item Uses efficient eigenvalue computation
}
}
\note{
Last change 15/03/2025
}
\keyword{internal}
