#' Generate Simulated Survival Accelerated Failure Time model Data with Two Views U and Z.
#'
#' @description
#' Function generate_data generates simulated survival data for Accelerated Failure Time (AFT) models with two views
#' here,  \code{U} and \code{Z}, that are correlated with unobserved latent variables (\code{L}).
#' The two observed views \code{U} and \code{Z} might represent two omics data matrices such as gene expression and methylation. The latent factor \code{L} can represent some unobserved factors that are associated with the survival times.
#' Some of the observed variables in  \code{U} and \code{Z} are correlated with the latent factors \code{L}.
#' The simulator allows to consider different proportion of right censoring.
#' The simulation scheme extends to the AFT models the simulation scheme proposed in Ding et al, 2022 for assessing the performance of cooperative regression,
#' and it is described in detail in Angelini et al 2025.
#' The function supports Weibull, Lognormal, and Loglogistic distributions for the error noise.
#'
#' @param model Character string specifying the noise distribution. In can assume one of the following values: "weibull", "lognormal", or "loglogistic".
#' @param n Integer; number of observations (it must be a positive integer number)
#' @param pu Integer; number of covariates in view \code{U} of dimension  \code{n x p_u}
#' @param pz Integer; number of covariates in view \code{Z} of dimension  \code{n x p_z}
#' @param tu Numeric; strength of view \code{U} in relationship with latent variables \code{L}.
#' @param tz Numeric; strength of view \code{Z} in relationship with latent variables \code{L}.
#' @param rate Numeric; minimum right censoring rate for the observed survival times.
#' @param sigma_true Numeric; true scale parameter for the AFT model
#' @param rl Integer; number of relevant latent covariates (default: 40).
#' @param rl_pos Integer; number of positively associated latent covariates (default: 20)
#' @param snr Numeric; signal-to-noise ratio (default: 0.8)
#'
#' @details
#' The data generation process follows these steps:
#' \itemize{
#'   \item Generates latent variables (L) from normal distribution
#'   \item Creates views U and Z with controlled relationship to \code{L}, to this purpose \code{tu} and \code{tz} represent the correlation among the rl relevant variables in the latent factor \code{L} and the first rl variables in both \code{U} and \code{Z}.
#'   \item Generates survival times based on specified distribution from the latent variables \code{L}.
#'   \item Applies random censoring to achieve desired censoring rate.
#' }
#' The function requires the flexsurv package for the loglogistic model.
#'
#' The function uses the following parameter ranges:
#' \itemize{
#'   \item Positive coefficients: Uniform(0.1, 1)
#'   \item Negative coefficients: Uniform(-1, -0.1)
#'   \item Censoring: Exponential distribution
#' }
#'
#' View strengths (tu, tz) typically use these combinations:
#' \itemize{
#'   \item (6,6): Strong and equal effect in both views U and Z.
#'   \item (6,4): Strong view in U, moderate effect in Z.
#'   \item (6,2): Strong view in U, weak effect in Z.
#'   \item (6,0): Strong view in U, No effect in Z.
#' }
#'
#' The choice of the SNR (default value snr=0.8) allow to generate different challenging scenarios.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item Y: Log-transformed censored survival times
#'   \item times_c: Original scale censored survival times
#'   \item delta: Censoring indicators (1 = event, 0 = censored)
#'   \item L: Matrix of latent variables
#'   \item beta_L: True coefficients for latent variables
#'   \item U: First view matrix
#'   \item Z: Second view matrix
#' }
#'
#'
#' @references
#' Angelini, De Canditiis, De Feis, Iuliano (in prep. 2025).
#'
#' Ding, D., Li, S., Narasimhan, B., Tibshirani, R. (2021) <doi:10.1073/pnas.2202113119>)
#'
#' @examples
#' \dontrun{
#' # Generate Weibull survival data
#' data <- generate_data(
#'   model = "weibull",
#'   n = 200,
#'   pu = 150,
#'   pz = 150,
#'   tu = 6,
#'   tz = 6,
#'   rate = 40,
#'   sigma_true = 0.5
#' )
#'
#' # Access components
#' Y <- data$Y  # log survival times
#' delta <- data$delta  # censoring indicators
#' U <- data$U  # first view
#' Z <- data$Z  # second view
#' }
#'
#' @import flexsurv
#'
#' @details
#' The function prints the actual percentage of non-censored observations and
#' a summary of the censored survival times. If the variance constraints cannot
#' be satisfied with the given parameters, it will adjust the view variances
#' and print a warning message.
#'
#' @note Last change 05/06/2025
#' @export

generate_data <- function(model,n,pu,pz,tu,tz,rate=40,sigma_true,rl=40,rl_pos=20,snr=0.8){

  library(flexsurv)

  # model="weibull";"lognormal";"loglogistic"
  # n=numbers of observations
  # pu=number of covariates on first view U
  # pz=number of covariates on first view Z
  # tu <- strenght of view U
  # tz <- strenght of view Z
  # (tu,tz)=(6,6);(6,4);(6,2);(6,0);
  # rate <- minimum censoring rate

  # parameters for the beta_L coefficients
  b1min <- 0.1
  b1max <- 1
  b2min <- -1
  b2max <- -0.1

  # if we want to mantain the scale with the same intensity
  # independently from the number of latent variables,
  # we have to rescale the variance according to rl
  # rl=40 #number of latent covariates
  # rl_pos=20 #number of positive associated latent covariates to the times, --> rl-rl_pos number of negatively associated latent covariates to the times

  sl  <- sqrt(1) #standard deviation for the latent

  # variance of the noise
  # (noise=standard minimum Gumbel/normal/logistic)
  if(model=="weibull"){
    var_standardnoise<-pi^2/6 ## variance of standard Weibull
  }
  if(model=="lognormal"){
    var_standardnoise<-1 ## variance of standard LogNormal
  }
  if(model=="loglogistic"){
    var_standardnoise<-pi^2/3 ## variance of standard LogLogistic
  }

  # standard deviation for U and Z:
  # i.e., the linear predictor f(X) accounted for
  # 90% of the variance of log(T)

  var_betaL_pos<-(b1max-b1min)^2/12 ## variance of uniform random variable in [b1min,b1max]
  mean_betaL_pos<-(b1max+b1min)/2  ## mean of uniform random variable in [b1min,b1max]
  var_betaL_neg<-(b2max-b2min)^2/12
  mean_betaL_neg<-(b2max+b2min)/2

  sum_betaL2<-rl_pos*(var_betaL_pos+mean_betaL_pos^2)+(rl-rl_pos)*(var_betaL_neg+mean_betaL_neg^2)

  # perc_var<-(1/0.9-1)

    perc_var <-(1/snr-1) ## to balance the Signal to noise

    if (sl^2*(tu+tz)^2*sum_betaL2*perc_var-sigma_true^2*var_standardnoise>0) {
      #   note when tu and tz are small, sluz cab be nagative, then we change sluz with a small value
      sluz <-sqrt((sl^2*(tu+tz)^2*sum_betaL2*perc_var-sigma_true^2*var_standardnoise)/((pu+pz)))
    } else {
      print(paste(" U and Z variance modified to guarantee positive contrains"))
      sluz<- 0.01 ## small variance
      # sluz <-sqrt((sl^2*(1)^2*sum_betaL2*perc_var-sigma_true^2*var_standardnoise)/((pu+pz)))
    }
  #
  U <- matrix(data = 0, nrow = n, ncol = pu, byrow = FALSE,dimnames = NULL)
  Z <- matrix(data = 0, nrow = n, ncol = pz, byrow = FALSE,dimnames = NULL)
  L <- matrix(data = 0, nrow = n, ncol = rl, byrow = FALSE,dimnames = NULL)

  for(j in 1:pu){
    U[,j] <- rnorm(n, mean = 0, sd = sluz) ## with noise sluz
  }
  for(j in 1:pz){
    Z[,j] <- rnorm(n, mean = 0, sd = sluz) ## with noise sluz
  }
  for(j in 1:rl){
    HVu <- NULL
    HVz <- NULL

    L[,j] <- rnorm(n, mean = 0, sd = sl)

    HVu <- U[,j]+tu *L[,j]
    U[,j] <- HVu

    HVz <- Z[,j]+tz *L[,j]
    Z[,j] <- HVz

  }

  beta_L <- rep(0, rl)
  beta_L[1:rl_pos] <- runif(rl_pos, min=b1min, max=b1max)
  beta_L[(rl_pos+1):rl] <- runif(rl-rl_pos, min=b2min, max=b2max)

  shape <- 1/sigma_true
  scale <- NULL
  # Generation of the survival time
  times <- vector()
  if(model=="weibull"){
    scale <- exp(L%*%beta_L)
    for(i in c(1:n)){
      times[i] <- rweibull(1, shape=shape, scale = scale[i])
    }
  }
  if(model=="lognormal"){
    scale <- L%*%beta_L
    for(i in c(1:n)){
      times[i] <- rlnorm(1, meanlog = scale[i], sdlog = sigma_true)
    }
  }
  if(model=="loglogistic"){
    scale <- exp(L%*%beta_L)
    for(i in c(1:n)){
      times[i] <- rllogis(1, shape = shape, scale = scale[i])
    }
  }
  ############################################################
  # generate censoring data from uniform with mean=theta=2/quantile(times,1-rate/100)
  # theta=2/quantile(times,1-rate/100)
  # censor <- runif(n,min=0,max=theta)
  #
  # generate censoring data from exponential with mean=1/theta=quantile(times,1-rate/100)
  theta=1/quantile(times,1-rate/100);
  delta <- rep(1,n)
  censor <-rexp(n,theta)
  times_c<-times
  for( i in c(1:n)){
    if (times[i]>censor[i]){
      delta[i]<-0
      times_c[i]<-censor[i]
    }
  }

  print(paste0("percentage not censored=",length(which(delta==1))/n))
  print(summary(times_c))
#  plot(times_c)

  Y<-log(times_c) ## times Y are already in log scale; times_c are in linear scale

  data <- list("Y"=Y,"times_c"=times_c,"delta"=delta,"L"=L,"beta_L"=beta_L,"U"=U,"Z"=Z)

}
