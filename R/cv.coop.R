#' K-fold Cross-validation for Cooperative AFT Model Estimation
#'
#' @description
#' Function cv.coop performs k-fold cross-validation for cooperative estimation in Accelerated Failure Time
#' (AFT) models. This function helps selecting the optimal tuning parameters (lambda) by
#' minimizing prediction error across folds. It is typically used by function \code{\link{aft_coop}} to select the optimal parameter.
#'
#' @param X Matrix of combined covariates of dimension \code{n x p_u+p_z} obtained as \code{X=[U,Z]}, when two views are provided (case=="coop"). Otherwise \code{X=[U]} or  \code{X=[Z]}, when case=="onlyU" or  case=="onlyZ", respectively.
#' @param Xtilde Matrix of cooperative adjusted covariates for cooperative estimation of dimension \code{n x p_u+p_z} obtained as \code{Xtilde=[U,-Z]}, when two views are provided (case=="coop");
#' @param Y Vector of log observed survival times or log censored times (response variable in log scale) of dimension \code{n}.
#' @param delta Vector of censoring indicators (1 for event, 0 for censored) of dimension \code{n}.
#' @param sigma Scale parameter for the AFT model.
#' @param lambda Vector of regularization parameters in descending order.
#' @param rho scalar cooperation parameter or vector of cooperation parameters, the latter option works when case=="coop".
#' @param nfolds Number of cross-validation folds (default: 5)
#' @param model Character string specifying the AFT model type. The following choices are available: "weibull";"lognormal";"loglogistic".
#' @param parallel Logical; whether to use parallel processing (default: TRUE)
#' @param ncore_max Maximum number of cores for parallel processing (default: 5)
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @details
#' The function implements k-fold cross-validation by:
#' \itemize{
#'   \item Randomly partitioning the data into k folds
#'   \item For each fold:
#'     \itemize{
#'       \item Fitting the model on training data (excluding the fold)
#'       \item Computing predictions on validation data (the fold)
#'       \item Calculating prediction error
#'     }
#'   \item Computing both minimum error lambda (lambda.min) and 1-standard-error rule
#'     lambda (lambda.1se)
#' }
#'
#' The function requires the cvTools package for fold generation.
#' The number of cores used for parallel processing is automatically determined
#' based on system availability, number of folds and user-specified maximum.
#' Parallel processing implementation differs between Windows (using makeCluster)
#' and Unix-like systems (using mclapply).
#'
#' @return
#' A list of class "cv.coop" containing:
#' \itemize{
#'   \item cv.err.linPred: Vector of cross-validation errors for each lambda
#'   \item cv.err.obj: Matrix of fold-specific errors (rows are folds, columns are lambdas)
#'   \item lambda_grid:  vector of values used in the grid of lambda for the cross validation
#'   \item lambda.min: Value of lambda that minimizes cross-validation error
#'   \item ind.lambda.min: Index of lambda.min in the lambda vector
#'   \item lambda.1se: Largest lambda within one standard error of minimum
#'   \item ind.lambda.1se: Index of lambda.1se in the lambda vector
#'   \item case.rho: chosed value of parameter rho
#' }
#'
#'
#' @import parallel
#' @importFrom cvTools cvFolds
#' @importFrom stats sd
#'
#' @seealso
#' \code{\link{aft_coop}} for the main estimation function
#'
#' @note Last change 25/03/2025
#'
#' @export
#'
#'
cv.coop <- function(X, Xtilde, Y, delta, sigma = NULL, lambda = NULL, rho = NULL,
                        nfolds = 5, model, parallel = TRUE, ncore_max=5,seed = 123) {

  set.seed(seed)
  library(parallel)

  n <- nrow(X)
  p <- ncol(X)

  # Create folds
  folds <- cvTools::cvFolds(n, K = nfolds, type = "random")
  nlambda <- length(lambda)  # Number of candidate tuning parameters to consider

  preds <- matrix(NA, nrow = n, ncol = nlambda)
  preds_coop <- matrix(NA, nrow = n, ncol = nlambda)
  cv.err.obj <- matrix(NA, nrow = nfolds, ncol = nlambda)

  out.fold <- vector("list", nfolds)  # To store fold CV results

  if (parallel == TRUE) {
    ncores <- min(c(detectCores() - 1, nfolds,ncore_max))

    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(ncores)
      clusterExport(cl, varlist = c("cv.aftcoop_fold","proxGD.coop","nll","hessian_coop", "gradient","prox.l1","folds", "X", "Xtilde", "Y", "delta",
                                    "sigma", "lambda", "rho", "model"), envir = environment())

      out.fold <- parLapply(cl, 1:nfolds, function(k) {
        cv.aftcoop_fold(folds, k, X, Xtilde, Y, delta, sigma = sigma, lambda = lambda, rho = rho, model = model)
      })

      stopCluster(cl)
    } else {
      out.fold <- mclapply(1:nfolds, function(k) {
        cv.aftcoop_fold(folds, k, X, Xtilde, Y, delta, sigma = sigma, lambda = lambda, rho = rho, model = model)
      }, mc.cores = ncores)
    }
  } else {
    for (k in seq_len(nfolds)) {
      out.fold[[k]] <- cv.aftcoop_fold(folds, k, X, Xtilde, Y, delta, sigma = sigma, lambda = lambda, rho = rho, model = model)
    }
  }

  # Combine results from parallel/sequential execution
  for (k in seq_len(nfolds)) {
    index.cv <- folds$subsets[folds$which != k]
    index.pred=setdiff(1:n ,index.cv)
    preds[index.pred, ] <- out.fold[[k]]$preds
    cv.err.obj[k, ] <- out.fold[[k]]$cv.err.obj
  }

  # Compute cross-validation error for each lambda more efficiently
  cv.err.linPred <- sapply(seq_len(nlambda), function(ll) {
    nll(Y, preds[, ll], delta, sigma, n, model)
  })

  # Find lambda.min
  ind_min <- which.min(cv.err.linPred)
  lambda.min <- lambda[ind_min]

  # Compute lambda.1se
  cv.err.linPred.se <- vapply(seq_len(nlambda), function(ll) {
    sd(cv.err.obj[, ll]) / sqrt(nfolds)
  }, numeric(1))

  lambda.1se <- max(lambda[cv.err.linPred < (cv.err.linPred[ind_min] + cv.err.linPred.se[ind_min])])
  ind_1se <- which(lambda==lambda.1se)

  # Output results
  out <- list("cv.err.linPred" = cv.err.linPred,
              "cv.err.obj" = cv.err.obj,
              "lambda_grid"=lambda,
              "lambda.min" = lambda.min,
              "ind.lambda.min" = ind_min,
              "lambda.1se" = lambda.1se,
              "ind.lambda.1se" = ind_1se,
              "case.rho"=rho)

  class(out) <- "cv.coop"
  return(out)
}

