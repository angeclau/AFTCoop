#' Compute Cross-validation for a Single Fold in AFT Cooperative Estimation
#'
#' @description
#' Internal function that performs cross-validation computations for a single fold
#' in the cooperative Accelerated Failure Time (AFT) model estimation. This function
#' handles data standardization, model fitting, and prediction error calculation
#' for one specific fold of k-fold cross-validation.
#'
#' @param folds List object containing fold information from cvTools::cvFolds
#' @param k Integer indicating the current fold number
#' @param X Matrix of combined covariates
#' @param Xtilde Matrix of cooperative term covariates for cooperative estimation
#' @param Y Vector of log observed survival times or log censored times (response variable in log scale)
#' @param delta Vector of censoring indicators (1 for event, 0 for censored)
#' @param sigma Scale parameter for the AFT model
#' @param lambda Vector of regularization parameters
#' @param rho Cooperation parameter
#' @param nfolds Number of cross-validation folds
#' @param model Character string specifying the AFT model type
#'
#' @details
#' The function performs the following steps for a single fold:
#' \itemize{
#'   \item Splits data into training and test sets based on fold indices
#'   \item Standardizes both X and Xtilde matrices using training set statistics
#'   \item Computes the cooperative term using standardized Xtilde
#'   \item For each lambda value:
#'     \itemize{
#'       \item Fits the model using proximal gradient descent
#'       \item Computes predictions on the test set
#'       \item Calculates prediction error
#'     }
#' }
#'
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item preds: Matrix of predictions for test set observations (rows) across
#'         different lambda values (columns)
#'   \item cv.err.obj: Vector of prediction errors for each lambda value
#' }
#'
#' @section Implementation Details:
#' \itemize{
#'   \item Uses efficient matrix operations (crossprod instead of t(X) %*% X)
#'   \item Implements warm starts for beta initialization across lambda values
#'   \item Handles standardization internally to ensure proper scaling
#' }
#'
#'
#' @seealso
#' \code{\link{cv.coop}} for the main cross-validation function
#' \code{\link{proxGD.coop}} for the underlying model fitting function
#'
#' @note Last change 03/04/2025
#' @keywords internal


cv.aftcoop_fold <- function(folds, k, X, Xtilde, Y, delta, sigma = NULL, lambda = NULL, rho = NULL,
                            nfolds = NULL, model) {

  # Get training indices
  index.cv <- folds$subsets[folds$which != k]
  index.test <- setdiff(seq_len(nrow(X)), index.cv)  # More efficient than `X[-index.cv,]`

  ntrain <- length(index.cv)
  ntest <- length(index.test)
  nlambda <- length(lambda)

  preds <- matrix(NA, nrow = ntest, ncol = nlambda)  # Preallocate

  # Standardization function
  standardize_with_mean_std <- function(data, mean_vals, sd_vals) {
    (data - matrix(mean_vals, nrow = nrow(data), ncol = length(mean_vals), byrow = TRUE)) /
      matrix(sd_vals, nrow = nrow(data), ncol = length(sd_vals), byrow = TRUE)
  }

  # Standardize training and test sets for X
  mX_train <- colMeans(X[index.cv, ])
  sdX_train <- apply(X[index.cv, ], 2, sd)
  X_train <-  standardize_with_mean_std(X[index.cv, ], mX_train, sdX_train)
  X_test <-  standardize_with_mean_std(X[index.test, ], mX_train, sdX_train) ##note use train mean and train std

  # Standardize training and test sets for Xtilde
  mXtilde_train <- colMeans(Xtilde[index.cv, ])
  sdXtilde_train <- apply(Xtilde[index.cv, ], 2, sd)
  Xtilde_train <- standardize_with_mean_std(Xtilde[index.cv, ], mXtilde_train, sdXtilde_train)
  Xtilde_test <- standardize_with_mean_std(Xtilde[index.test, ], mXtilde_train, sdXtilde_train) #note use train mean and train std


  # Precompute cooperative term
  coop_train <- crossprod(Xtilde_train)  # More efficient than t(Xtilde_train) %*% Xtilde_train

  # Extract response variables
  Y_train <- Y[index.cv]
  delta_train <- delta[index.cv]
  Y_test <- Y[index.test]
  delta_test <- delta[index.test]

  # Initialize beta0
  beta0 <- rep(0, ncol(X_train))
 cv.err.obj<-vector()
  # Cross-validated predictions loop
  for (ll in seq_len(nlambda)) {
    cv.getPath <- proxGD.coop(X_train, Y_train, delta_train, coop_train, beta0, sigma, rho, lambda[ll], model)
    beta0 <- cv.getPath$beta  # Warm start update
    preds[, ll] <- X_test %*% beta0  # Linear prediction
     cv.err.obj[ll] <- nll(Y_test, preds[, ll], delta_test, sigma, length(Y_test), model)
  }

  return(list(preds = preds, cv.err.obj = cv.err.obj))
}
