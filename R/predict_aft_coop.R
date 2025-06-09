#' Evaluation of the linear predictor score for the atf_coop Survival Model
#'
#' This function computes the linear predictor for the cooperative AFT survival model,
#' where the linear predictor is defined as \eqn{X \hat{\beta}} with \eqn{X} being the regression matrix and  \eqn{\hat{\beta}} the vector of estimated parameters obtained from the \code{\link{aft_coop}} function.
#' The linear predictor score can be used as Prognostic index to stratify patients in risk classes or to compute C-Index.
#'
#' @param U A matrix of feature variables (e.g., gene expression data). Each row corresponds to a sample, and each column corresponds to a feature (gene). It should have dimensions \code{n x pu}.
#' @param Z A matrix of additional feature variables (e.g., methylation data). Each row corresponds to a sample, and each column corresponds to a feature (attribute). It should have dimensions \code{n x pz}.
#' @param mU a vector of parameters for centering the U matrix. (default mU=colMeans(U), but when predicting on the test set the value mus be adjusted to those of the training set.)
#' @param mZ a vector of parameters for centering the Z matrix. (default mZ=colMeans(Z),but when predicting on the test set the value mus be adjusted to those of the training set.)
#' @param beta A vector of estimated coefficients (of length (pu), (pz) or (pu + pz), depending on the case) obtained from obtained from the \code{\link{aft_coop}} function.
#' @param case A character string specifying the case for prediction. Options are:
#' \itemize{
#' \item{"onlyU"}{Predict using only the  \code{U} matrix.}
#' \item{"onlyZ"}{Predict using only the \code{Z} matrix.}
#' \item{"coop"}{Predict using both the  \code{U} and \code{Z}  matrices combined.}
#'    }
#'
#' @details
#' The depending on the estimation strategy of the regression coefficients, the prediction can be made based on:
#' \itemize{
#'   \item "onlyU": When only the \code{U} matrix (e.g., gene expression data) was used to estimate \eqn{\hat{\beta}}.
#'   \item "onlyZ": When only the \code{Z} matrix (e.g., methylation data) was used to estimate \eqn{\hat{\beta}}.
#'   \item "coop": When  the two matrices \code{U} and \code{Z} matrices were used to estimate \eqn{\hat{\beta}}.
#' }
#'
#' Moreover, when applied to the test set the regression matrices are centered according to the mean of the training set.
#' This is achieved using the parameter mU and mZ.
#'
#'
#' @return A vector of predicted linear scores.
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
#'     rho_values<- c(1,0,0.25,0.5,0.75) #vector parameteters for rho values
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
#' ## linear prediction using the estimated coefficients corresponding to the first value of rho_values
#' lin.pred<- predict_aft_coop(U=U, Z=Z, beta=beta_est_coop[,1], case = "coop")
#' print(lin.pred)  # Predicted values based on combined data
#' #' }
#'
#' @note Last change 11/03/2025
#' @export
#'
predict_aft_coop <- function(U, Z, mU=NULL, mZ=NULL,beta, case=c("onlyU", "onlyZ","coop")){

  if(case=="onlyU"){
    if (is.null(mU)){
      mU <- colMeans(U)
    }
    n <- dim(U)[1]
    UC <-  U - tcrossprod(rep(1, n), mU) # center the data
    pred_vals <- UC %*% beta
  }

  if(case=="onlyZ"){
    if (is.null(mZ)){
      mZ <- colMeans(Z)
    }
    n <- dim(Z)[1]
    ZC <-  Z - tcrossprod(rep(1, n), mZ) # center the data
    pred_vals <- ZC %*% beta
  }

  if(case=="coop"){
    if (is.null(mU) & is.null(mZ)){
      mU <- colMeans(U)
      mZ <- colMeans(Z)
    }
    UZ <- cbind(U,Z)
    n <- dim(UZ)[1]
    mUZ <- c(mU, mZ)
    UZC <- UZ - tcrossprod(rep(1, n), mUZ)
    pred_vals <- UZC %*% beta
  }

  return(pred_vals=pred_vals)

}
