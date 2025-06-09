#' Plot AFT Cooperative Cross-Validation Results
#'
#' @description
#' Internal function that creates a ggplot visualization of cross-validation results
#' from an accelerated failure time (AFT) cooperative regularization model.
#'
#' @param cv.out List containing cross-validation results obtained by \code{\link{cv.coop}} , must include components:
#'   \itemize{
#'     \item cv.err.obj - Matrix of CV errors for each fold and lambda value
#'     \item cv.err.linPred - Vector of mean CV errors for linear predictor
#'     \item  lambda_grid Numeric vector of lambda values used in the regularization path.
#'     \item lambda.1se - Lambda value within 1 standard error of minimum CV error
#'     \item lambda.min - Lambda value with minimum CV error
#'     \item case.rho - Rho value for the current case
#'   }
#' @param nfolds Integer specifying the number of cross-validation folds used.
#' @param case Character or numeric identifier for the current case being analyzed.
#' @param addbar Logical to add error bar to the crossvalidation curve (Default FALSE)
#'
#' @details
#' The function calculates standard errors for cross-validation error estimates
#' and creates a visualization with error bars. It marks the lambda.min and lambda.1se
#' values with vertical dashed lines. The object cv.err.obj is obtained using \code{\link{cv.coop}}
#' @return A ggplot object visualizing cross-validation results
#' @import ggplot2
#'
#'
#' @note Last change 15/04/2025
#'
#'
#' @keywords internal


plot_aft_coop.cv<-function(cv.out,nfolds,case,addbar=F){
  library(ggplot2)
 # print("Plotting CV function")
  cv.err.linPred.se <- apply(cv.out$cv.err.obj,2,sd)/sqrt(nfolds)
  cvup <-   cv.out$cv.err.linPred+cv.err.linPred.se
  cvlo <-  cv.out$cv.err.linPred-cv.err.linPred.se
  df = data.frame(lambda = log(cv.out$lambda_grid), LPCVE = cv.out$cv.err.linPred, LPCVEhi=cvup, LPCVElow=cvlo)
  plot <- ggplot(df,aes(x=lambda,y=LPCVE)) +
    geom_point(col="#f05454") +
    xlab(expression(log(lambda))) +
    geom_vline(xintercept=log(cv.out$lambda.1se),linetype="dotted",col="blue") +
    geom_vline(xintercept=log(cv.out$lambda.min),linetype="dashed",col="orange") +
    labs(y = "Linear predictor CV error ") + ggtitle(paste0("CV cross validation for case ",case, " with rho= ",cv.out$case.rho))+
    theme_bw()
  if (addbar){
    plot<-plot+geom_errorbar(aes(ymin = LPCVElow,ymax=LPCVEhi),col="#30475e")
  }

  print(plot)
}
