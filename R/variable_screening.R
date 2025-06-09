#' Variable Screening method
#'
#' @description
#' This function selects a subset of relevant features from a set of potential features based on their association with patient survival,
#' thus allowing a preliminary reduction in the number of variables before applying a penalized regression approach.
#' The function uses the following ranking methods:
#' \itemize{
#' \item \code{absmg}: The marginal utilities ordered from largest to smallest.
#' \item \code{mg}: The regression coefficients ordered from largest to smallest, with the top and bottom \code{perc} selected.
#' \item \code{mgpadj}: The marginal utilities with p-values adjusted to be less than \code{thresh}.
#' }
#' @param X Matrix of covariates correspondig to the view of dimension \code{n x p_{view}}, where \code{n} is the number of samples and \code{p_{view}} is the number of variables in the view.
#' @param y Vector of observed survival times or censored times of dimension \code{n}, where \code{n} is the number of samples.
#' @param delta Vector of censoring indicators (1 for event, 0 for censored) of dimension  \code{n}, where \code{n} is the number of samples.
#' @param family Character string indicating the type of survival model to fit. Options are \code{cox} for Cox Proportional Hazards model and \code{AFT} for Accelerated Failure Time model.
#' @param model Character string specifying the specific AFT model to use.  It can assume one of the following values "weibull";"lognormal";"loglogistic". Ignored if family = "cox".
#' @param rank Character string specifying the ranking method.  It can assume one of the following values:"absmg", "mg" or "mgpadj"
#' @param topn An integer specifying the number of top features to return (default 500). Used when rank = "absmg".
#' @param perc An integer specifying the percentage of top and bottom features to return (default  5). Used when rank = "mg".
#' @param thresh A numeric value representing the threshold for adjusted p-values. Only features with adjusted p-values smaller than \code{thresh} will be considered. Used when rank = "mgpadj".
#'
#' @details
#' The function performs variable screening based on the marginal ranking of the features in matrix \code{X}, thus reducing the features' dimensionality while retaining relevant ones.
#' :
#' \itemize{
#' \item \code{absmg}: It uses the absolute marginal utility function and selects the \code{topn} features as screened variables.
#' \item \code{mg}:  It uses the marginal utility function and selects as screened variables given percentage \code{perc} of the top or bottom features.
#' \item \code{mgpadj}: It uses the p-value adjusted according to the marginal utility function to select as screened variables those with adjusted p-values less than  \code{thresh}.
#' }
#' The adjustment method includes the Benjamini and Hochberg (1995) (\code{BH} or its alias \code{fdr}) correction.
#'
#' @return A vector of indices corresponding to the \code{topn} ranked features when rank = "absmg",
#' the top and bottom features when rank = "mg", and the features with adjusted p-values smaller than \code{thresh} when rank = "mgpadj".
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' model <- "weibull"
#' topn <- 500
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
#' y <- data$times_c    # vector of survival times or censored times
#' delta <- data$delta  # vector of censoring indicators
#' U <- data$U  # matrix corresponding to the first view
#' colnames(U) <- paste0("feature_",  1:ncol(U))
#' Z <- data$Z  # matrix corresponding to the second view
#' colnames(Z) <- paste0("feature_",  1:ncol(Z))
#'
#' index_screen_U <- variable_screening(U, y, delta, family = "AFT", model = model, rank = "absmg", topn = topn)
#' index_screen_Z <- variable_screening(Z, y, delta, family = "AFT", model = model, rank = "absmg", topn = topn)
#'}
#'
#'
#' @seealso
#' \code{\link{marginal_ranking}} for the evaluation of the marginal ranking of the individual features.
#'
#' @note Last change 08/04/2025
#' @export

variable_screening <- function(X, y, delta,family = c("cox", "AFT"), model = c("weibull", "loglogistic", "lognormal"), rank = c("absmg", "mg", "mgpadj"), topn = 500, perc = 5, thresh = 0.5){

  # Ranking by marginal utility
  screen_feature <- marginal_ranking(X, y, delta, family = family, model = model, rank = rank)

  if(rank=="absmg"){
    index <- screen_feature$ranking[1:topn]
  }

  if(rank=="mg"){
    top_ud<-round(perc*length(screen_feature$ranking))
    topup <- head(screen_feature$ranking, top_ud)
    topdown <- tail(screen_feature$ranking,top_ud)
    index <- c(topup,topdown)
  }

  if(rank=="mgpadj"){
    index <- screen_feature$ranking[which(screen_feature$mgpadj < thresh)]
  }

  return(index=index)

}
