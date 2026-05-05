#' Variable Screening method
#'
#' @description
#' This function selects a subset of relevant features from a set of potential features based on their association with patient survival,
#' allowing a preliminary reduction in the number of variables before applying a penalized regression approach.
#'
#' The function supports the following ranking methods:
#' \itemize{
#' \item \code{absmg}: Marginal utilities ordered from largest to smallest (absolute value).
#' \item \code{mg}: Regression coefficients ordered from largest to smallest; both top and bottom \code{perc} are selected.
#' \item \code{mgpadj}: Marginal utilities filtered by adjusted p-values less than \code{thresh}.
#' \item \code{mgpval}: Variables ranked using \eqn{-\log_{10}(p\text{-value}) \times |\beta|}, combining significance and effect size.
#' }
#'
#' @param X Matrix of covariates of dimension \code{n x p_view}, where \code{n} is the number of samples and \code{p_view} is the number of variables.
#' @param y Vector of observed survival or censoring times of length \code{n}.
#' @param delta Vector of censoring indicators (1 = event, 0 = censored) of length \code{n}.
#' @param family Character string specifying the survival model. Options are \code{cox} or \code{AFT}.
#' @param model Character string specifying the AFT model. Options are \code{weibull}, \code{lognormal}, or \code{loglogistic}. Ignored if \code{family = cox}.
#' @param rank Character string specifying the ranking method. Options are \code{absmg}, \code{mg}, \code{mgpadj}, or \code{mgpval}.
#' @param topn Integer specifying the number of top features to return (default is 500). Used when \code{rank = absmg} or \code{rank = mgpval}.
#' @param perc Numeric value in \code{(0,1)} specifying the proportion of top and bottom features to select (default is 0.1). Used when \code{rank = mg}.
#' @param thresh Numeric threshold for adjusted p-values (default is 0.5). Used when \code{rank = mgpadj}.
#'
#' @details
#' The function performs variable screening based on marginal ranking of features in \code{X}, reducing dimensionality while retaining relevant variables.
#' \itemize{
#' \item \code{absmg}: Selects the top \code{topn} features based on absolute marginal utility.
#' \item \code{mg}: Selects both the top and bottom \code{perc} proportion of features based on marginal effects.
#' \item \code{mgpadj}: Selects features with adjusted p-values below \code{thresh}, using the Benjamini-Hochberg correction (\code{BH} or \code{fdr}).
#' \item \code{mgpval}: Selects the top \code{topn} features ranked by a combined score of p-values and effect sizes.
#' }
#'
#' @return
#' A vector of indices corresponding to the selected features:
#' \itemize{
#' \item Top \code{topn} features for \code{absmg} and \code{mgpval}.
#' \item Top and bottom features for \code{mg}.
#' \item Features with adjusted p-values below \code{thresh} for \code{mgpadj}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' model <- "weibull"
#' topn <- 500
#'
#' data <- generate_data(
#'   model = model,
#'   n = 200,
#'   pu = 150,
#'   pz = 150,
#'   tu = 6,
#'   tz = 6,
#'   rate = 40,
#'   sigma_true = 0.5
#' )
#'
#' y <- data$times_c
#' delta <- data$delta
#'
#' U <- data$U
#' colnames(U) <- paste0("feature_", 1:ncol(U))
#'
#' Z <- data$Z
#' colnames(Z) <- paste0("feature_", 1:ncol(Z))
#'
#' index_screen_U <- variable_screening(
#'   U, y, delta,
#'   family = "AFT",
#'   model = model,
#'   rank = "absmg",
#'   topn = topn
#' )
#'
#' index_screen_Z <- variable_screening(
#'   Z, y, delta,
#'   family = "AFT",
#'   model = model,
#'   rank = "absmg",
#'   topn = topn
#' )
#' }
#'
#' @seealso
#' \code{\link{marginal_ranking}} for computing marginal feature rankings.
#'
#' @note Last updated: 2026-04-22
#'
#' @export

variable_screening <- function(X, y, delta,family = c("cox", "AFT"), model = c("weibull", "loglogistic", "lognormal"), rank = c("absmg", "mg", "mgpadj", "mgpval"), topn = 500, perc = 0.1, thresh = 0.5){

  # Ranking by marginal utility
  screen_feature <- marginal_ranking(X, y, delta, family = family, model = model, rank = rank)

  if(rank=="absmg"){
    index <- screen_feature$ranking[1:topn]
  }

  if (rank == "mg") {
    r <- screen_feature$ranking
    top_ud <- round(perc * length(r))
    topup <- head(seq_along(r), top_ud)
    topdown <- tail(seq_along(r), top_ud)
    index <- c(topup, topdown)
  }

  if(rank == "mgpadj"){
    index <- screen_feature$ranking[which(screen_feature$mgpadj < thresh)]
  }

  if (rank == "mgpval") {
    index <- screen_feature$ranking[1:topn]
  }

  return(index = index)

}
