% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_ranking.R
\name{marginal_ranking}
\alias{marginal_ranking}
\title{Variable screening based on ranking by marginal utility}
\usage{
marginal_ranking(
  X,
  y,
  delta,
  family = c("cox", "AFT"),
  model = c("weibull", "exponential", "lognormal", "loglogistic"),
  rank = c("absmg", "mg", "mgpadj")
)
}
\arguments{
\item{X}{Matrix of covariates correspondig to the view of dimension \code{n x p_{view}}, where \code{n} is the number of samples and \code{p_{view}} is the number of variables in the view.}

\item{y}{Vector of observed survival times in the time-domain or censored times of dimension \code{n}, where \code{n} is the number of samples.}

\item{delta}{Vector of censoring indicators (1 for event, 0 for censored) of dimension  \code{n}, where \code{n} is the number of samples.}

\item{family}{Character string indicating the type of survival model to fit. Options are \code{cox} for Cox Proportional Hazards model and \code{AFT} for Accelerated Failure Time model.}

\item{model}{Character string specifying the specific AFT model to use.  It can assume one of the following values: "weibull";"lognormal";"loglogistic". Ignored if \code{family = "cox"}.}

\item{rank}{Character string specifying the ranking method. It can assume one of the following values: "absmg", "mg" and "mgpadj".}
}
\value{
A data frame with the ranked features. The columns include:
\itemize{
\item ranking:  The rank of the feature based on the selected method.
\item symbol: The name of the feature.
\item rank: The ranked marginal utility values using one of several methods: \code{absmg}, \code{mg}, \code{mgpadj}.
}
}
\description{
This function performs marginal ranking of features based on their association with the survival outcome
(see, for example, Fan and Lv (2008) and Fan, Feng and Wu (2010)).
Specifically, it computes the marginal coefficient for each feature in the dataset and s
orts the covariates in decreasing order according to the selected ranking method:
\code{absmg}, \code{mg}, or \code{mgpadj}.
}
\details{
This function implements marginal ranking for the features of a matrix \code{X} based on their association with survival outcomes in \code{y,delta} ,
using either the Cox Proportional Hazards or Accelerated Failure Time (AFT) models.
The functions assume that the feature in \code{X}  were already preprocessed to remove irrelevant information and normalized, if necessary.
Such preprocessing depends on the type of variables encoded in matrix \code{X}.
Then, it uses the \code{margcoef} function to compute the marginal coefficients for each feature,
employing several ranking methods (namely,\code{absmg}, \code{mg}, or \code{mgpadj}) based on marginal utility or adjusted p-values.
The \code{margcoef} function, in turn, uses the \code{mg} function to calculate the regression coefficients
and adjusted p-values, depending on the specified \code{family} (\code{Cox} or \code{AFT}).
The ranking of the features can be done according to one of the following criteria:
\itemize{
\item \code{absmg}: The absolute marginal utility function.
\item \code{mg}: The marginal utility function.
\item \code{mgpadj}: The p-value adjusted according with the marginal utility function.
}
When choosing \code{mgpadj}, the adjustment method includes the Benjamini and Hochberg (1995) (\code{BH} or its alias \code{fdr}) correction.
}
\note{
Last change 08/04/2025
}
\examples{
\dontrun{
# Generate example data
set.seed(123)
model <- "weibull"
data <- generate_data(
  model = model,
    n = 200,
    pu = 150,
    pz = 150,
    tu = 6,
    tz = 6,
    rate = 40,
    sigma_true = 0.5
    )

y <- data$times_c    # vector of survival times or censored times
delta <- data$delta  # vector of censoring indicators
U <- data$U  # matrix corresponding to the first view
colnames(U) <- paste0("feature_",  1:ncol(U))
Z <- data$Z  # matrix corresponding to the second view
colnames(Z) <- paste0("feature_",  1:ncol(Z))

screening_U <- marginal_ranking(U, y, delta, family = "AFT", model = model, rank = "absmg")
print(screening_U)
screening_Z <- marginal_ranking(U, y, delta, family = "AFT", model = model, rank = "absmg")
print(screening_U)
}

}
\references{
Benjamini, Y., and Hochberg, Y. (1995). "Controlling the false discovery rate: a practical and powerful approach to multiple testing",
Journal of the Royal Statistical Society Series B, 57, 289–300.
DOI: 10.1111/j.2517-6161.1995.tb02031.x

Fan, J., Feng, Y., and Wu, Y. (2010). "High-dimensional Variable Selection for Cox Proportional Hazards Model".
IMS Collections, 6, 70-86.
DOI: 10.1214/10-IMSCOLL606

Fan, J., and Lv, J. (2008). "Sure independence screening for ultrahigh dimensional feature space",
Journal of the Royal Statistical Society Series B: Statistical Methodology, 70(5), 849-911.
DOI: 10.1111/j.1467-9868.2008.00674.x
}
