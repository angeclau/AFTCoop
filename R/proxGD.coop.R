#' Proximal Gradient Descent for Cooperative AFT Model Estimation
#'
#' @description
#' proxGD.coop implements proximal gradient descent optimization for cooperative estimation in
#' Accelerated Failure Time (AFT) models. The function minimizes an objective function
#' that combines negative log-likelihood, L1 regularization, and a cooperative term.
#'
#' @param X Matrix of standardized covariates
#' @param Y Vector of log observed survival times or log censored times (response variable in log scale)
#' @param delta Vector of censoring indicators (1 for event, 0 for censored)
#' @param coop Cooperative matrix (typically t(Xtilde) %*% Xtilde)
#' @param beta0 Initial vector of coefficients
#' @param sigma Scale parameter for the AFT model
#' @param rho Cooperation parameter controlling L1 vs cooperative trade-off (
#' @param lambda Regularization parameter
#' @param model Character string specifying the AFT model type
#' @param niter Maximum number of iterations (default: 1000)
#'
#' @details
#' The function minimizes the objective:
#' \deqn{f(\beta) = -\log L(\beta) + \lambda\rho\|\beta\|_1 + \lambda(1-\rho)\beta^T C\beta}
#' where:
#' \itemize{
#'   \item \eqn{-\log L(\beta)} is the negative log-likelihood
#'   \item \eqn{\|\beta\|_1} is the L1 norm (lasso penalty)
#'   \item \eqn{\beta^T C\beta} is the cooperative term
#' }
#'
#' The algorithm uses:
#' \itemize{
#'   \item Adaptive line search with backtracking
#'   \item Lipschitz constant adaptation
#'   \item Early stopping based on relative change in objective
#'   \item Memory-efficient vector operations
#' }
#'
#' Convergence is determined by relative change in objective value falling below
#' \code{conv = 0.001}.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item beta: Vector of optimized coefficients
#'   \item objective: Final value of the objective function
#' }
#'
#' @section Implementation Notes:
#' \itemize{
#'   \item Uses pre-allocation for memory efficiency
#'   \item Implements line search with sufficient decrease condition
#'   \item Precomputes lambda terms to avoid repeated multiplication
#'   \item Uses vectorized operations where possible
#' }
#'
#'
#' @seealso
#' \code{\link{prox.l1}} for the proximal operator of L1 norm
#' \code{\link{nll}} for negative log-likelihood computation
#' \code{\link{gradient}} for gradient computation
#' \code{\link{hessian_coop}} for Hessian and Lipschitz constant computation
#'
#'
#' @importFrom stats coef
#' @note Last change 07/04/2025
#' @keywords internal

proxGD.coop <- function(X, Y, delta, coop, beta0, sigma = NULL, rho = NULL,
                       lambda = NULL, model, niter = 1000) {
  # Pre-compute constants
  n <- length(Y)
  conv <- 0.001
  value <- 2
  p <- ncol(X)

  # Pre-allocate vectors instead of matrix for memory efficiency
  beta_current <- beta0
  beta_next <- numeric(p)

  # Initialize objective tracking and use drop() to ensure proper vector dimensionality
  eta <- as.vector(X %*% beta_current)
  part0 <- drop(t(beta_current) %*% coop %*% beta_current)
  negloglike_0 <- nll(Y, eta, delta, sigma, n, model)
  obj_current <- negloglike_0 + lambda * rho * sum(abs(beta_current)) +
                 lambda * (1 - rho) * part0

  # Initial Lipschitz constant
  L <- hessian_coop(X, Y, eta, delta, sigma, rho, lambda, coop, n, model)$L
  fun_old <- negloglike_0 + lambda * (1 - rho) * part0

  # Precompute lambda terms
  lambda_rho <- lambda * rho
  lambda_1_rho <- lambda * (1 - rho)

  # Main optimization loop
  for (k in 1:niter) {
    # Compute gradient
    gradient_beta <- gradient(X, Y, eta, delta, sigma, n, model)$grad.beta
    coop_beta <- coop %*% beta_current

    # Line search loop
    repeat {
      # Proximal gradient step
      u_beta <- beta_current - (1/L) * (gradient_beta + 2 * lambda_1_rho * coop_beta)
      beta_next <- prox.l1(u_beta, lambda_rho/L)

      # Update function values
      eta_new <- as.vector(X %*% beta_next)
      part1 <- drop(t(beta_next) %*% coop %*% beta_next)
      negloglike_new <- nll(Y, eta_new, delta, sigma, n, model)
      fun_new <- negloglike_new + lambda_1_rho * part1

      # Check sufficient decrease condition
      beta_diff <- beta_next - beta_current
      if (fun_new <= fun_old +
          drop(t(gradient_beta + 2 * lambda_1_rho * coop_beta) %*% beta_diff) +
          (L/2) * sum(beta_diff^2)) {
        break
      }

      L <- value * L

    }

    # Update objective
    obj_new <- negloglike_new + lambda_rho * sum(abs(beta_next)) +
               lambda_1_rho * part1

    # Check convergence
    rel_change <- abs(obj_new - obj_current)/abs(obj_current)
    if (rel_change <= conv) {
      break
    }

    # Update states for next iteration
    beta_current <- beta_next
    eta <- eta_new
    obj_current <- obj_new
    fun_old <- fun_new
  }
  list(beta = beta_next, objective = obj_new)

}
