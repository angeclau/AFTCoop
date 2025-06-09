#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Negative log-likelihood function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Calculate Negative Log-Likelihood for AFT Models
#'
#' @description
#' Computes the negative log-likelihood for Accelerated Failure Time (AFT) models,
#' supporting Weibull, lognormal, and loglogistic distributions.
#'
#' @param Y Vector of observed times (response variable)
#' @param eta Vector of linear predictors
#' @param delta Vector of censoring indicators (1 for event, 0 for censored)
#' @param sigma Scale parameter
#' @param n Sample size
#' @param model Character string specifying the distribution ("weibull", "lognormal", or "loglogistic")
#'
#' @return
#' Normalized negative log-likelihood value (i.e., divided by sample size)
#'
#'
#' @seealso \code{\link{ll}} for the positive log-likelihood version
#'
#' @note Last change 15/03/2025
#' @keywords internal
#'
nll <- function(Y,eta,delta,sigma,n,model){

  res <- (Y-eta)/sigma

  loss<-switch(model,
               "weibull" = as.numeric(crossprod(delta, (-log(sigma) + res)) - sum(exp(res))),
               "lognormal" = as.numeric(-1/2*(log(2*pi*sigma^2))*sum(delta)-1/2*crossprod(delta,res^2)+crossprod((1-delta),log(1-pnorm(res)))) ,
               "loglogistic" = as.numeric(crossprod(delta,(-log(sigma) + res))-crossprod((1+delta),log1p(exp(res)))),
               stop("Invalid model specified")  # Handles unexpected model names
  )
  loglike <- -loss/n
  return(loglike)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Gradient of negative log-likelihood function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Calculate Gradient of Negative Log-Likelihood for AFT Models
#'
#' @description
#' Computes the gradient of the negative log-likelihood function for Accelerated
#' Failure Time (AFT) models with respect to the coefficients.
#'
#' @param X Matrix of covariates
#' @inheritParams nll
#'
#' @details
#' Implements efficient matrix operations and pre-computations for better performance.
#' Gradient calculations are specific to each supported distribution type.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item grad.beta: Vector of gradients with respect to coefficients
#' }
#'
#' @note Last change 24/03/2025
#' @keywords internal

gradient <- function(X, Y, eta, delta, sigma, n, model) {

  res <- (Y - eta) / sigma
  exp_res <- exp(res)  # Precompute exp(res) once to avoid redundant calculations

  a <- switch(model,
              "weibull" = exp_res - delta,
              "lognormal" = (delta * res) + ((1 - delta) * dnorm(res) / (1 - pnorm(res))),
              "loglogistic" = (exp_res / (1 + exp_res))* (1 + delta) - delta,
              stop("Invalid model specified")  # Handles unexpected model names
  )

  grad.beta <- - (1 / (n * sigma)) * crossprod(X, a)  # More efficient than t(X) %*% a

  return(list(grad.beta = grad.beta))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Hessian plus Lipschitz constant of negative log-likelihood function plus coop penalty
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' Calculate Hessian and Lipschitz Constant for AFT Models
#'
#' @description
#' Computes the Hessian matrix and its maximum eigenvalue (Lipschitz constant)
#' for the negative log-likelihood plus cooperative penalty.
#'
#' @param X Matrix of covariates
#' @param rho Cooperation parameter
#' @param lambda Regularization parameter
#' @param coop Cooperative matrix
#' @inheritParams nll
#' @details
#' The function:
#' \itemize{
#'   \item Computes second derivatives specific to each distribution
#'   \item Includes cooperative penalty term in Hessian
#'   \item Uses efficient eigenvalue computation
#' }
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item L: Lipschitz constant (maximum eigenvalue of Hessian)
#' }
#' @note Last change 15/03/2025
#' @keywords internal
hessian_coop <- function(X, Y, eta, delta, sigma, rho, lambda, coop, n, model) {
  res <- (Y - eta) / sigma  # Compute residuals once

  E <- switch(model,
              "weibull" = diag(exp(res)),
              "lognormal" = {
                dens <- dnorm(res)
                surv <- 1 - pnorm(res)
                diag(delta + (1 - delta) * (dens * (-res * surv + dens)) / surv^2)
              },
              "loglogistic" = {
                exp_res <- exp(res)
                diag((1 + delta) * exp_res / ((1 + exp_res)^2))
              },
              stop("Invalid model specified")  # Handles unexpected model names
  )

  H <- (1 / n) * crossprod(X, E %*% X) + lambda * (1 - rho) * coop
  max.eig <- max(eigen(H, only.values = TRUE)$values)  # Get max eigenvalue directly

  return(list(L = max.eig))  # Return Lipschitz constant
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Proximal operator: lasso
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#' L1 Proximal Operator
#'
#' @description
#' Computes the proximal operator for L1 regularization (soft thresholding).
#'
#' @param u Numeric vector of input values
#' @param mu Threshold parameter (typically lambda * alpha / L)
#'
#' @details
#' Implements the soft thresholding operator:
#' \deqn{prox_{\mu\|\cdot\|_1}(u) = sign(u)\max(|u| - \mu, 0)}
#'
#' @return
#' Vector of proximal operator values
#'
#' @note Last change 11/03/2025
#' @keywords internal

prox.l1 <- function(u, mu) { # mu = lambda*alpha/L
  uhat <- abs(u) - mu
  ubind <- cbind(rep(0, length(u)), uhat)
  prox <- sign(u) * apply(ubind, 1, max)
  return(prox)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Matrix standartization
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##' Standardize Matrix Columns
#'
#' @description
#' Standardizes a matrix by centering and scaling each column to have zero mean
#' and unit variance.
#'
#' @param mat Input matrix to be standardized
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item std: Standardized matrix
#'   \item means: Vector of column means
#'   \item sds: Vector of column standard deviations
#' }
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), 20, 5)
#' std_result <- standardize(X)
#' # Access standardized matrix
#' X_std <- std_result$std
#' # Original scale restoration
#' X_restored <- scale(X_std, center = -std_result$means, scale = 1/std_result$sds)
#' }
#'
#' @note Last change 11/03/2025
#' @keywords internal
#'
standardize <- function(mat) {
  scaled_mat <- scale(mat)  # Standardizes columns
  return(list(
    std = scaled_mat,
    means = attr(scaled_mat, "scaled:center"),
    sds = attr(scaled_mat, "scaled:scale")
  ))
}

