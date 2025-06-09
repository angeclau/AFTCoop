#' Cooperative Estimation with AFT Survival Models
#'
#' @description
#' The function aft_coop implements a cooperative estimation method for Accelerated Failure Time (AFT) models
#' allowing to integrate two matrices \code{U} and \code{Z} of high-dimensional covariates.
#' The two matrices  \code{U} and \code{Z} can represent two omics views, such as gene expression and methylation.
#' When used with two matrices  \code{U} and \code{Z}, the function estimates the coefficients \eqn{\beta_u} and \eqn{\beta_z}
#' associated with the two views by solving the minimization problem using the proximal gradient descend algorithm.
#' In such a case, the minimization problem consists of three terms: the negative log-likelihood, the agreement term among the two views,
#' and the Lasso penalization for each view.
#' However, the function works also with a single matrix  \code{U} or \code{Z}.
#' In such cases, it performs negative log-likelihood minimization with the Lasso penalty.
#' In both cases, cross-validation estimates the regularization parameter.
#' The function supports parallel processing and cross-validation for optimal parameter selection.
#'
#' @param U Matrix of covariates correspondig to the first view (e.g., gene expression) of dimension \code{n x p_u}, where \code{n} is the number of samples and \code{p_u} is the number of variables in the first view.
#' @param Z Matrix of covariates correspondig to the second view (e.g., gene methylation) of dimension \code{n x p_z},where \code{n} is the number of samples and \code{p_z} is the number of variables in the second view.
#' @param Y Vector of log observed survival times or log censored times (i.e., the response variables are in log scale) of dimension \code{n}, where \code{n} is the number of samples.
#' @param delta Vector of censoring indicators (1 for event, 0 for censored) of dimension  \code{n}, where \code{n} is the number of samples.
#' @param sigma Scale parameter for the AFT model (it must be a positive real value).
#' @param nfolds Number of folds to use for the cross-validation (default: 5)
#' @param model Character string specifying the specific AFT model to use.  It can assume one of the following values "weibull";"lognormal";"loglogistic".
#' The noise distribution in the log time models  \eqn{Y=X \beta+\sigma \epsilon} can be Extreme Value Gumbel distribution that corresponds to the "weibull" model, normal distribution that correspond to the  "lognormal" model, and logistic distribution that correspond to the "loglogistic" model.
#' @param case Character string specifying the data model. It can take the following values: "onlyU","onlyZ","coop". The choice "onlyU" uses only the first view of covariates,
#'             The choice "onlyZ" uses only the second view of covariates, The choice "coop" uses both views as covariates.
#' @param rho_values Vector of \code{nrhos} cooperation parameters. The  \code{nrhos} parameters are relevant only when case="coop" and allow different tuning of the coperative and Lasso penalization terms. When case="onlyU" or case="onlyZ", the parameter is not relevant and it is set to 1.
#' @param lam_min Logical; if TRUE uses lambda.min, if FALSE uses lambda.1se
#' @param nlambda Number of lambda values in regularization path (default: 100)
#' @param lambda.ratio.min Minimum ratio of smallest to largest lambda (default: 0.01)
#' @param parallel Logical; enable parallel processing (default: TRUE)
#' @param ncore_max_rho Maximum number of cores for parallel processing for splitting over rho_values (default: 4)
#' @param ncore_max_cv Maximum number of cores for parallel processing for splitting over Cross-validation (default: 5)
#' @param seed Random seed for reproducibility (default: 123)
#' @param iplot Logical; enable the plot of the CV function (default: TRUE)
#' @param addbar Logical; when iplot is true, add the error bar to the plot (default: FALSE)
#'
#' @details
#' When providing two matrices \code{U} and \code{Z}, the function aft_coop performs the following steps:
#' \itemize{
#'   \item Build combined matrix \code{X=[U;Z]} and \code{Xtilde=[U;-Z]}, given the two matrices \code{U} and \code{Z}
#'   \item Standardizes \code{U} and \code{Z}  matrices.
#'   \item Generates lambda grid and initialize the search using gradient at initial estimate \eqn{\beta = 0}.
#'   \item For each rho value:
#'     \itemize{
#'       \item Performs K-fold cross-validation to select optimal lambda
#'       \item Computes cooperative matrix
#'       \item Applies proximal gradient descent
#'       \item Rescales coefficients back to original scale
#'     }
#' }
#' The function uses standardized variables internally but returns coefficients
#' in the original scale. Cross-validation is used to select optimal lambda values.
#' The user can choose the optimal value for lambda as lambda.min (i.e., the values of lambda for which the CV is minimised)
#' or lambda.1se (i.e, the largests value of lambda within 1 standard error form the the minumum).
#' The latter choice provides a sparsest model.
#' The choice of the type of optimal value to use is done setting the logical lam_min parameter.
#' However, we switch between  lambda.min  and lambda.1se is the correspondoing values is at the boarder of the grid interval.
#' Parallel processing implementation differs between Windows (using makeCluster)
#' and Unix-like systems (using mclapply).
#'
#' @return
#' A matrix of dimension \code{(p_u+p_z)x nrhos},where:
#' \itemize{
#'   \item Rows correspond to the estimates regression coefficients \eqn{\hat{\beta}=[\hat{\beta}_u,\hat{\beta}_z]} for each view.
#'   \item Columns correspond to the estimates for different rho values
#'   \item Values are coefficient \eqn{\hat{\beta}} in original scale.
#'   The first set of coefficients corresponds \eqn{\hat{\beta}_u} to those associated with the \code{U} view, the second set \eqn{\hat{\beta_z}} to those associated to the view  \code{Z}.
#' }
#'
#'
#' @references
#' Angelini, De Canditiis, De Feis, Iuliano (in prep. 2025).
#'
#' Ding, Li, Narasimhan, Tibshirani, "Cooperative learning for multiview analysis", PNAS (2022)
#' DOI: 10.1073/pnas.220211311
#' Available at: \url{https://www.pnas.org/doi/abs/10.1073/pnas.2202113119}
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' model="weibull"
#' data <- generate_data(
#'   model = model,
#'     n = 200,
#'     pu = 150,
#'     pz = 150,
#'     tu = 6,
#'     tz = 6,
#'     rate = 40,
#'     sigma_true = 0.5
#'     )
#'
#'     Y <- data$Y  # vector of log survival times or censored times
#'     delta <- data$delta  # vector of censoring indicators
#'     U <- data$U  # matrix corresponding to the first view
#'     Z <- data$Z  # matrix corresponding to the second view
#'     fit_survreg <- survreg(Surv(data$times_c, delta)~1, dist=model,scale=0)
#'     sigma.est <- exp(fit_survreg$icoef[2])
#'     rho_values<- c(1,0.25,0.5,0.75) #vector parameteters for rho values
#' # Using aft_coop with two views
#'     beta_est_coop<-aft_coop(
#'     U=U,
#'     Z=Z,
#'     Y=Y,
#'     delta=delta,
#'     sigma=sigma.est,
#'     nfolds=5,
#'     model=model,
#'     case="coop",
#'     rho_values=  rho_values,
#'     lam_min=F,
#'     parallel=T,
#'     ncore_max=5,
#'     seed=123)
#'
#'    ## estimate of beta  (the first pu components refers to the view U; The second pz component to view Z  )
#'    beta_est_coop
#'
#' # Using aft_coop with only the first view U
#'     beta_est_onlyU<-aft_coop(
#'     U=U,
#'     Z=Z,
#'     Y=Y,
#'     delta=delta,
#'     sigma=sigma.est,
#'     nfolds=5,
#'     model=model,
#'     case="onlyU",
#'     rho_values=  1,
#'     lam_min=F,
#'     parallel=T,
#'     ncore_max=5,
#'     seed=123)
#'
#'
#' ## Using aft_coop with only the second view
#'  beta_est_onlyZ<-aft_coop(
#'     U=U,
#'     Z=Z,
#'     Y=Y,
#'     delta=delta,
#'     sigma=sigma.est,
#'     nfolds=5,
#'     model=model,
#'     case="onlyZ",
#'     rho_values=  1,
#'     lam_min=F,
#'     parallel=T,
#'     ncore_max=5,
#'     seed=123)
#'
#'     beta_est_onlyZ
#' }
#'
#'
#' @import parallel ggplot2
#' @importFrom stats rnorm rbinom
#'
#'
#' @seealso
#' \code{\link{cv.coop}} for underlying implementation details
#'
#' @note Last change 15/04/2025
#' @export


aft_coop <- function(U, Z, Y, delta, sigma, nfolds=5, model, case, rho_values,
                         lam_min, nlambda = 100, lambda.ratio.min = 0.01, parallel = TRUE,
                         ncore_max_rho = 4,ncore_max_cv=5, seed = 123, iplot=T,addbar=F) {


  set.seed(seed)
  library(parallel)
  library(ggplot2)

  # Check input types and values
  if (!is.matrix(U) || !is.numeric(U)) stop("U must be a numeric matrix.")
  if (!is.matrix(Z) || !is.numeric(Z)) stop("Z must be a numeric matrix.")
  if (!is.vector(Y) || !is.numeric(Y) || length(Y) != nrow(U)) stop("Y must be a numeric vector of length n (number of samples).")
  if (!is.vector(delta) || !is.numeric(delta) || length(delta) != nrow(U)) stop("delta must be a numeric vector of length n (number of samples).")
  if (!is.numeric(sigma) || sigma <= 0) stop("sigma must be a positive real number.")
  if (!is.numeric(nfolds) || nfolds <= 0 || nfolds != round(nfolds)) stop("nfolds must be a positive integer.")
  if (!is.character(model) || !(model %in% c("weibull", "lognormal", "loglogistic"))) stop("model must be one of 'weibull', 'lognormal', or 'loglogistic'.")
  if (!is.character(case) || !(case %in% c("onlyU", "onlyZ", "coop"))) stop("case must be one of 'onlyU', 'onlyZ', or 'coop'.")
  if (!is.numeric(rho_values) || length(rho_values) == 0) stop("rho_values must be a numeric vector with at least one element. The elements must be non negatives.")
  if (!is.logical(lam_min)) stop("lam_min must be a logical value.")
  if (!is.numeric(nlambda) || nlambda <= 0 || nlambda != round(nlambda)) stop("nlambda must be a positive integer.")
  if (!is.numeric(lambda.ratio.min) || lambda.ratio.min <= 0 || lambda.ratio.min >= 1) stop("lambda.ratio.min must be a number between 0 and 1.")
  if (!is.logical(parallel)) stop("parallel must be a logical value.")
  if (!is.numeric(ncore_max_rho) || ncore_max_rho <= 0 || ncore_max_rho != round(ncore_max_rho)) stop("ncore_max_rho must be a positive integer.")
  if (!is.numeric(ncore_max_cv) || ncore_max_cv <= 0 || ncore_max_cv != round(ncore_max_cv)) stop("ncore_max_cv must be a positive integer.")
  if (!is.numeric(seed) || seed != round(seed)) stop("seed must be an integer.")


  # Set Dimensions
  pu <- ncol(U)
  pz <- ncol(Z)
  n <- length(Y)

  n_rho_values <- length(rho_values)

  # Standardize
  U_std_res <- standardize(U)
  Z_std_res <- standardize(Z)

  # Select case
  if (case == "onlyU") {
    p <- pu
    X<-U
    Xtilde<-U
    X_std <- U_std_res$std
    Xtilde_std <- X_std
  } else if (case == "onlyZ") {
    p <- pz
    X<-Z
    Xtilde<-Z
    X_std <- Z_std_res$std
    Xtilde_std <- X_std
  } else {
    p <- pu + pz
    X<-cbind(U,Z)
    Xtilde<-cbind(U,-Z)
    X_std <- cbind(U_std_res$std, Z_std_res$std)
    Xtilde_std <- cbind(U_std_res$std, -Z_std_res$std)
  }

  # Initialize and Compute lambda grid
  beta0 <- rep(0, p)
  eta <- X_std %*% beta0
  grad_beta <- gradient(X_std, Y, eta, delta, sigma, n, model)$grad.beta

  # Initialize cooperative matrix penalty
  coop <- crossprod(Xtilde_std, Xtilde_std)  # t(Xtilde_std) %*% Xtilde_std

  print("----------------------------------")
  print(paste("Running AFTcoop algorithm for case:", case))
  print("----------------------------------")

  # Output matrix
  beta_est <- matrix(0, nrow = p, ncol = n_rho_values)

  # Parallel processing
  process_rho <- function(rho) {
      lambda_max <- max(abs(grad_beta))/rho
      lambda_min <- lambda.ratio.min * lambda_max
      lambda_grid <- exp(seq(log(lambda_max),log(lambda_min),length.out=nlambda))
      cv_out <- cv.coop(X, Xtilde, Y, delta, sigma, lambda = lambda_grid, rho, nfolds =nfolds, model=model,
                          parallel=parallel, ncore_max = ncore_max_cv, seed = seed + round(100 * rho))

## we avoid to take lambda equal to lambda_min or lambda_max
if (lam_min) {
   if (cv_out$lambda.min> lambda_min ){
     lambda <- cv_out$lambda.min
   } else {lambda <-cv_out$lambda.1se }
 } else {
   if (cv_out$lambda.1se<lambda_max){
     lambda <-cv_out$lambda.1se
  }    else {lambda <-cv_out$lambda.min}
}

    out <- proxGD.coop(X_std, Y, delta, coop, beta0, sigma, rho, lambda, model)
    res_beta<-out$beta / if (case == "onlyU") U_std_res$sds else if (case == "onlyZ") Z_std_res$sds else c(U_std_res$sds, Z_std_res$sds)
    res<-list(beta=res_beta,cv_out= cv_out)
    return(res)
  }

  if (parallel && length(rho_values) > 1) {
    ncores <- min(detectCores() - 1, length(rho_values), ncore_max_rho)

      all_results <- if (.Platform$OS.type == "windows") {
      cl <- makeCluster(ncores)
      sdU=U_std_res$sds
      sdZ=Z_std_res$sds
      clusterExport(cl, varlist = c("cv.coop", "proxGD.coop","cv.aftcoop_fold","nll","gradient","hessian_coop","prox.l1",
                                    "X", "Xtilde", "Y", "delta", "sigma", "grad_beta", "model","nfolds", "parallel","seed", "lambda.ratio.min","nlambda",
                                    "case", "sdU", "sdZ","coop","ncore_max_cv"), envir = environment())
      results <- parLapply(cl, rho_values, process_rho)
      stopCluster(cl)
      results
    } else {
      mclapply(rho_values, process_rho, mc.cores = ncores)
    }

    for (jj in 1:length(all_results)){
      beta_est[,jj]<-all_results[[jj]]$beta
      cv.out<-all_results[[jj]]$cv_out
      if (iplot==T){
        plot_aft_coop.cv(cv.out,nfolds,case,addbar)
      }
    }
  } else {
    for (i in seq_along(rho_values)) {
      all_results<- process_rho(rho_values[i])
      beta_est[, i] <-all_results$beta
      cv.out<-all_results$cv_out
      if (iplot==T){
        plot_aft_coop.cv(cv.out,nfolds,case,addbar)
      }
    }
  }
  return(beta_est)
}
