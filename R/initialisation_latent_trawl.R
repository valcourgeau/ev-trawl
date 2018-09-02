#' Zeta function from auto-correlation expansion for small lags.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter. Must be positive.
#' @param kappa Exceedance probability parameter. Must be positive.
#' @examples Zeta(alpha = 5, beta = 4, kappa = 3)
#' @export
Zeta <- function(alpha, beta, kappa){
  res.zeta <- beta^alpha / ((beta+kappa)^{alpha-1})
  res.zeta <- res.zeta * hypergeo::hypergeo(A = alpha-1, B = alpha-1, C = alpha, z = - beta/(beta+kappa))
  return(Re(res.zeta))
}

#' Computes first-order estimate for \code{rho}.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter. Must be positive.
#' @param kappa Exceedance probability parameter. Must be positive.
#' @param cluster.size Lag at which ACF is negligeable i.e. cluster size.
#' @param data Exceedane timeseries to be used.
#'
#' @examples
#' # TODO ADD data
#' GetEstimateRho(alpha = 5, beta = 2, kappa = 3, cluster.size = 8, data = NA)
#'
#' @export
GetEstimateRho <- function(alpha, beta, kappa, cluster.size, data){
  requireNamespace("stats", quietly = TRUE)
  if(alpha < 0){
    data[data>0] <- vapply(data[data>0], function(x){TrfInverseG(x, alpha = alpha,
                                                               beta = beta, kappa = kappa,
                                                               offset_scale = 2,
                                                               offset_shape = 1+kappa)}, FUN.VALUE = 1.0)
    alpha <- 4
    beta <- 1
  }

  d_plus <- beta^2 / ((alpha-2)*(alpha-1))*(1+2*kappa/beta)^{2-alpha}
  d_plus <- d_plus * (log(1+2*kappa/beta) + (2*alpha-3)/((alpha-2)*(alpha-1)))

  d_times <- 2*beta/((alpha-2)*(alpha-1))*(beta*(1+2*kappa/beta)^{2-alpha}*log(1+kappa/beta)+Zeta(alpha = alpha-1, beta = beta, kappa = kappa))
  return( abs(stats::var(data) / (cluster.size * alpha * (d_plus-d_times))))
}

#' Computes initial guess for Univariate Latent-Trawl model.
#' Uses \code{fExtreme} package method \code{gpdFit} to fit MLE.
#'
#' @param data Exceedane timeseries to be used.
#' @param cluster.size Lag at which ACF is negligeable.
#'
#' @examples
#' # TODO ADD data
#' data(hourlyhourly_bloomsbury_air_pollution_2000_2017)
#' GenerateParameters(data = hourly_bloomsbury_air_pollution_2000_2017$O3[1:1000], cluster.size = 8)
#'
#' @export
GenerateParameters <- function(data, cluster.size){
  requireNamespace("fExtremes", quietly = TRUE)
  requireNamespace("stats", quietly = TRUE)

  params_to_work_with <- rep(0, 4)
  fit_marginal <-  fExtremes::gpdFit(data[data > 0], u= 0)@fit$fit$par
  p_nz <- length(which(data > 0))/length(data)

  params_to_work_with <- rep(0, 4)
  params_to_work_with[1] <- 1/fit_marginal[1]
  params_to_work_with[2] <- abs(fit_marginal[2]*params_to_work_with[1])

  params_to_work_with[4] <- abs(params_to_work_with[2]*(p_nz^{-1/params_to_work_with[1]}-1) / (p_nz^{-1/params_to_work_with[1]}))
  params_to_work_with[2] <- abs(params_to_work_with[2] - params_to_work_with[4])

  params_to_work_with[3] <- GetEstimateRho(alpha = params_to_work_with[1],
                                             beta =  params_to_work_with[2],
                                             kappa = params_to_work_with[4],
                                             cluster.size,
                                             data)

  if(is.na(params_to_work_with[3])){
    # if rho is not a valid number, uses random initialisation.
    warning("rho initial value is not valid.")
    params_to_work_with[3] <- stats::runif(min = 1e-6, max = 1, n = 1)
  }
  return(params_to_work_with)
}
