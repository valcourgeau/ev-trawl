#' Wrapper to compute total area under exponential trawl function.
#'
#' @param rho Exponential trawl parameter (should be non-negative).
#'
#' @return Total area under exponential trawl function (1/rho).
#'
#' @examples
#' ComputeAExp(1) # should be 1.0
#' ComputeAExp(0.2) # should be 5.0.
ComputeAExp <- function(rho){
  if(rho < 0) stop ('rho should be non-negative.')
  return(1 / rho)
}

#' Wrapper to compute difference area between \code{t1} and \code{t2}
#' under exponential trawl function.
#'
#' @param rho Exponential trawl parameter (should be non-negative).
#' @param t1 First timestamp.
#' @param t2 Second timestamp.
#' @return Difference area between \code{t1} and \code{t2}
#' under exponential trawl function.
#'
#' @examples
#' ComputeB3Exp(1, t1 = 3, t2 = 5)
#' ComputeB3Exp(0.2, t1 = 7, t2 = 3)
ComputeB3Exp <- function(rho, t1, t2){
  if(rho < 0) stop('rho should be non-negative.')

  if(t2 > t1){
    return((1 - exp(rho * (t1-t2)))/rho)
  }else{
    if(t2 == t1){
      return(0.0)
    }else{
      return((1 - exp(rho * (t2 - t1)))/rho)
    }
  }
}

#' Wrapper to compute difference area between \code{t2} and \code{t1}
#' under exponential trawl function.
#'
#' @param rho Exponential trawl parameter (should be non-negative).
#' @param t1 First timestamp.
#' @param t2 Second timestamp.
#' @return Difference area between \code{t1} and \code{t2}
#' under exponential trawl function. Relies on ComputeB3Exp.
#'
#' @examples
#' ComputeB1Exp(1, t1 = 3, t2 = 5)
#' ComputeB1Exp(0.2, t1 = 7, t2 = 3)
ComputeB1Exp <- function(rho, t1, t2){
  if(rho < 0) stop('rho should be non-negative.')
  # Compute A_t_1 \ A_t_2
  return(ComputeB3Exp(rho, t2, t1))
}

#' Wrapper to compute intersection area between \code{t1} and \code{t2}
#' under exponential trawl function.
#'
#' @param rho Exponential trawl parameter (should be non-negative).
#' @param t1 First timestamp.
#' @param t2 Second timestamp.
#' @return Intersection area between \code{t1} and \code{t2}
#' under exponential trawl function.
#'
#' @examples
#' ComputeBInterExp(1, t1 = 3, t2 = 5)
#' ComputeBInterExp(0.2, t1 = 7, t2 = 3)
ComputeBInterExp <- function(rho, t1, t2){
  if(rho < 0) stop('rho should be non-negative.')

  if(t1 > t2){
    return(exp(rho * (t2-t1))/rho)
  }else{
    if(t1 == t2){
        return(1/rho)
    }else{
      return(exp(rho * (t1-t2))/rho)
    }
  }
}

#' Computes inverse of g from the marginal transform method (vectorised version).
#' That is from GPD(offset_shape, offset_scale) to GPD(alpha, beta+kappa).
#'
#' @param z Data distributed as GPD(offset_shape, offset_scale).
#' @param alpha Shape parameter of output data.
#' @param beta Part of scale parameter of output data (beta + kappa).
#' Should be positive.
#' @param kappa Part of scale parameter of output data.
#' Should be positive.
#' @param offset_scale Scale parameter of input data.
#' Should be positive.
#' @param offset_shape Shape parameter of input data.
#'
#' @return GPD(alpha, beta+kappa)-distributed data from
#' GPD(offset_shape, offset_scale)-distributed input data.
#'
#' @example
#' TrfInverseG(c(0.1, 0.2), 2, 3, 2, 4, 1)
TrfInverseG <- function(z, alpha, beta, kappa, offset_scale, offset_shape){
  # From GPD(alpha, beta+kappa) to GPD(offset_shape, offset_scale)
  if(beta <= 0) stop('beta should be positive.')
  if(kappa <= 0) stop('kappa should be positive.')
  if(offset_scale <= 0) stop('offset_scale should be positive.')

  res <- (offset_scale)*((1+sign(alpha)*z/(beta+kappa))^{alpha/offset_shape}-1)

  return(res)
}

#' Computes g from the marginal transform method (vectorised version).
#' That is from GPD(alpha, beta+kappa) to GPD(offset_shape, offset_scale).
#'
#' @param z Data distributed as GPD(offset_shape, offset_scale).
#' @param alpha Shape parameter of output data.
#' @param beta Part of scale parameter of output data (beta + kappa).
#' Should be positive.
#' @param kappa Part of scale parameter of output data.
#' Should be positive.
#' @param offset_scale Scale parameter of input data.
#' Should be positive.
#' @param offset_shape Shape parameter of input data.
#'
#' @return GPD(offset_shape, offset_scale)-distributed data from
#' GPD(alpha, beta+kappa)-distributed data input data.
#'
#' @example
#' TrfG(c(0.1, 0.2), 2, 3, 2, 4, 1)
TrfG <- function(x, alpha, beta, kappa, offset_scale, offset_shape){
  # From GPD(offset_shape, offset_scale) to GPD(alpha, beta+kappa)
  if(beta <= 0) stop('beta should be positive.')
  if(kappa <= 0) stop('kappa should be positive.')
  if(offset_scale <= 0) stop('offset_scale should be positive.')

  res <- sign(alpha)*(beta+kappa)*((1+sign(offset_shape)*x/(offset_scale))^{offset_shape/alpha}-1)
  return(res)
}

#' Computes predefined transformed scale parameter.
#'
#' @param alpha Shape parameter of output data.
#' @param beta Part of scale parameter of output data (beta + kappa).
#' Should be positive.
#' @param kappa Part of scale parameter of output data.
#' Should be positive.
#' @param offset_shape Shape parameter of input data.
#'
#' @return Transformed scale parameter given original data parameters.
#'
#' @example
#' TrfFindOffsetScale(c(0.1, 0.2), 2, 3, 2, 4, 1)
TrfFindOffsetScale <- function(alpha, beta, kappa, offset_shape){
  if(beta <= 0) stop('beta should be positive.')
  if(kappa <= 0) stop('kappa should be positive.')

  return(1+kappa)
}

#' Computes the probability density function of Generalised Pareto Distribution at x,
#' with shape parameter alpha and scale parameter beta.
#'
#' @param x value at which the pdf is evaluated.
#' @param alpha Shape parameter.
#' @param beta Scale parameter, should be positive.
#'
#' @return GPD pdf function evaluated at x with shape and scale parameters respectively alpha and beta.
#' @example dlgpd(2.34, alpha = 1, beta = 2)
dlgpd <- function(x, alpha, beta){
  if(beta <= 0) stop('beta should be positive.')
  return(abs(alpha)/beta*max(0,(1+sign(alpha)*x/beta))^{-alpha-1.0})
}

#' Computes the cumulative distribution function of Generalised Pareto Distribution at x,
#' with shape parameter alpha and scale parameter beta.
#'
#' @param x value at which the pdf is evaluated.
#' @param alpha Shape parameter.
#' @param beta Scale parameter, should be positive.
#'
#' @return GPD CDF function evaluated at x with shape and scale parameters respectively alpha and beta.
#' @example plgpd(2.34, alpha = 1, beta = 2)
plgpd <- function(x, alpha, beta, lower.tail=F){
  res <- 1-(1+x/beta)^{-alpha}
  if(lower.tail){
    res <- 1-res
  }
  return(res)
}

#' Computes determinant of Jacobian of Marginal Transform (MT) method.
#'
#' @param z GPD(offset_shape, offset_scale)-distributed data.
#' @param alpha Shape parameter of output data.
#' @param beta Part of scale parameter of output data (beta + kappa).
#' Should be positive.
#' @param kappa Part of scale parameter of output data.
#' Should be positive.
#' @param offset_scale Scale parameter of input data.
#' Should be positive.
#' @param offset_shape Shape parameter of input data.
#' @return Transformed scale parameter given original data parameters.
#'
#' @example
#' TrfJacobian(c(0.1, 0.2), 2, 3, 2, 1, 4)
TrfJacobian <- function(z, alpha, beta, kappa, offset_scale, offset_shape){
  # TODO check whether it is numerically stable by division of pdfs
  inv_g_z <- TrfInverseG(z = z, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
  res <- dlgpd(x = z, alpha = alpha, beta = beta+kappa) / dlgpd(x = inv_g_z, alpha = offset_shape, beta = offset_scale)
  return(res)
}

# Example
n_sample <- 1000
alpha <- 6
beta <- 12
kappa <- 4
proba_trf_k <- 1/(1+kappa)
proba_trf_a_b <- (1+1/beta)^{-alpha}

proba_no_trf <- (1+kappa/beta)^{-alpha}
library(fExtremes)
library(evir)

# TODO rewrite example

# Case 0-0

#' Computes partial part of latent trawl pairwise likelihood with \code{(x,y) = (0,0)}.
#'
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#'
#' @return Partial part of latent trawl pairwise likelihood with \code{(x,y) = (0,0)}.
#' @example PairwiseZeroZero1(0.3, 2, 3)
PairwiseZeroZero1 <- function(alpha, beta, kappa){
  return(-2 * (1 + kappa / beta)^{-alpha})
}

#' Computes second part of latent trawl pairwise likelihood with \code{(x,y) = (0,0)}.
#'
#' @param t1 First timestamp.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B1 Difference area between times \code{t1} and \code{t2}.
#' @param B2 Intersection area between times \code{t1} and \code{t2}.
#' @param B3 Difference area between times \code{t2} and \code{t1}.
#'
#' @return Second part of latent trawl pairwise likelihood with \code{(x,y) = (0,0)}.
#' @example PairwiseZeroZero2(t1=1, t2=4, 0.3, 2, 3, 0.2, 1.2, 3, 1.2)
PairwiseZeroZero2 <- function(t1, t2, alpha, beta, kappa, rho, B1, B2, B3){
  # TODO Remove dependency on rho
  # TODO Remove either B1 or B3
  temp = (1 + 2 * kappa / beta)^{-alpha * rho * B2}
  return((1 + kappa / beta)^{-alpha * rho * (B1 + B3)} * temp)
}

#' Computes latent trawl pairwise likelihood with \code{(x,y) = (0,0)} for Exponential
#' Trawl funciton.
#'
#' @param t1 First timestamp.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be non-negative.
#'
#' @return Second part of latent trawl pairwise likelihood with \code{(x,y) = (0,0)}.
#' @example PairwiseZeroZeroExp(t1=1, t2=4, 0.3, 2, 3, 0.2)
PairwiseZeroZeroExp <- function(t1, t2, alpha, beta, kappa, rho){
  if(rho < 0) stop('rho should be non-negative')
  B1 <- ComputeB1Exp(rho, t1, t2)
  B2 <- ComputeBInterExp(rho, t1, t2)
  B3 <- ComputeB3Exp(rho, t1, t2)

  temp <- PairwiseZeroZero1(alpha, beta, kappa)
  temp <- temp + PairwiseZeroZero2(t1, t2, alpha, beta, kappa, rho, B1, B2, B3)
  return(1 + temp)
}

# Example
t1 <- 0.0
t2 <- 1.0
alpha <- -1.0
beta <- 10
rho <- 1.0
kappa <- 1.0
B1 <- ComputeB1Exp(rho, t1, t2)
B2 <- ComputeBInterExp(rho, t1, t2)
B3 <- ComputeB3Exp(rho, t1, t2)

PairwiseZeroZeroExp(t1 = t1, t2 = t2,
                alpha = alpha, beta = beta,
                kappa = kappa, rho = rho)
answer <- 1 - 2 * (1 + 0.1) + (1 + 0.1)^{B1 + B3}*(1 + 0.2)^{B2}

alpha <- 1.0
PairwiseZeroZeroExp(t1 = t1, t2 = t2,
                alpha = alpha, beta = beta,
                kappa = kappa, rho = rho)
answer <- 1 - 2 * (1 + 0.1)^{-1} + (1 + 0.1)^{-B1 - B3}*(1 + 0.2)^{-B2}
answer

# Case 1-0

# TODO Remove dependency on rho

#' Computes first part of latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param trawlA Total trawl set measure / area under trawl function.
#'
#' @return first part of latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#' @example PairwiseOneZero1(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3)
PairwiseOneZero1 <- function(t1, x1, t2, alpha, beta, kappa, rho, trawlA){
  return(alpha * rho * trawlA / beta * (1 + (kappa + x1) / beta)^{-alpha * rho * trawlA - 1})
}

#' Computes partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B1 Difference area between \code{t2} and \code{t1} (in this order).
#'
#' @return Partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#' @example PairwiseOneZero21(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3)
PairwiseOneZero21 <- function(t1, x1, t2, alpha, beta, kappa, rho, B1){
  return( - alpha * rho / beta * (1 + (kappa + x1) / beta)^{-alpha * rho * B1 - 1})
}

#' Computes partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B2 Intersection area between \code{t1} and \code{t2}.
#'
#' @return Partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#' @example PairwiseOneZero22(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3)
PairwiseOneZero22<- function(t1, x1, t2, alpha, beta, kappa, rho, B2){
  return((1 + (2*kappa + x1) / beta)^{-alpha * rho * B2 - 1})
}

#' Computes partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B3 Difference area between \code{t1} and \code{t2} (in this order).
#'
#' @return Partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#' @example PairwiseOneZero23(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3)
PairwiseOneZero23 <- function(t1, x1, t2, alpha, beta, kappa, rho, B3){
  return((1 + kappa / beta)^{-alpha * rho * B3})
}

#' Computes partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param trawlA Total trawl set measure / area under trawl function.
#' @param B3 Difference area between \code{t1} and \code{t2} (in this order).
#'
#' @return Partial part of second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#' @example PairwiseOneZero24(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3, 3, 2)
PairwiseOneZero24 <- function(t1, x1, t2, alpha, beta, kappa, rho, trawlA, B1){
  # TODO remove dependcy on rho
  return(trawlA * (1 + (kappa + x1) / beta) + B1 * kappa / beta)
}

#' Computes second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param trawlA Total trawl set measure / area under trawl function.
#' @param B1 Difference area between \code{t1} and \code{t2} (in this order).
#' @param B2 intersection area between \code{t1} and \code{t2} (in this order).
#' @param B3 Difference area between \code{t2} and \code{t1} (in this order).
#'
#' @return Second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}.
#' @example PairwiseOneZero2(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3, 3, 2, 1, 2)
PairwiseOneZero2 <- function(t1, x1, t2, alpha, beta, kappa, rho, trawlA, B1, B2, B3){
  temp <- PairwiseOneZero21(t1, x1, t2, alpha, beta, kappa, rho, B1)
  temp <- temp * PairwiseOneZero22(t1, x1, t2, alpha, beta, kappa, rho, B2)
  temp <- temp * PairwiseOneZero23(t1, x1, t2, alpha, beta, kappa, rho, B3)
  temp <- temp * PairwiseOneZero24(t1, x1, t2, alpha, beta, kappa, rho, trawlA, B1)

  return(temp)
}

#' Computes term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}
#' with exponential trawl function.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param transformation Boolean to use the Marginal Transform (MT) method.
#' @param n_moments Number of moments achieved by transformed GPD marginals, if used.
#'
#' @return Second term in latent trawl pairwise likelihood with \code{(x,0)} where \code{x > 0}
#' with exponential trawl function.
#' @example PairwiseOneZeroExp(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3, 3, 2)
PairwiseOneZeroExp <- function(t1, x1, t2, alpha, beta, kappa, rho, transformation=F, n_moments=4){
  # Marginal Transformation
  if(transformation){
    offset_shape <- n_moments + 1
    offset_scale <- TrfFindOffsetScale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
    inv_x <- TrfInverseG(x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    jacobian <- TrfJacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    new_x <- inv_x
  }else{
    new_x <- x1
  }

  trawlA <- ComputeAExp(rho)
  B1 <- ComputeB1Exp(rho, t1, t2)
  B2 <- ComputeBInterExp(rho, t1, t2)
  B3 <- ComputeB3Exp(rho, t1, t2)

  if(transformation){
    temp <- PairwiseOneZero1(t1, new_x, t2, alpha=offset_shape, beta=1, kappa, rho, trawlA)
    temp <- temp + PairwiseOneZero2(t1, new_x, t2, alpha=offset_shape, beta=1, kappa, rho, trawlA, B1, B2, B3)
    temp <- temp * jacobian
  }else{
    temp <- PairwiseOneZero1(t1, new_x, t2, alpha, beta, kappa, rho, trawlA)
    temp <- temp + PairwiseOneZero2(t1, new_x, t2, alpha, beta, kappa, rho, trawlA, B1, B2, B3)
  }

  if(temp == 0.0 || is.na(temp) || is.nan(temp)){
    temp <- 1.0
  }

  return(temp)
}

# Example

# Case 1-1

#' Computes frst part of first term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param x2 Positive value corresponding to \code{t2}.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B1 Difference area between \code{t1} and \code{t2} (in this order).
#'
#' @return Frst part of first term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#' @example PairwiseOneOne11(t1=1, x1=0.5, t2=4, x2=1.5, 0.3, 2, 3, 0.2, 3, 2)
PairwiseOneOne11 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1){
  return(alpha^2 * rho^2 / beta^2 * (1+(kappa+x1)/beta)^{-alpha*rho*B1-1})
}

#' Computes second part of first term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param x2 Positive value corresponding to \code{t2}.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B2 Intersection area between \code{t1} and \code{t2} (in this order).
#'
#' @return Second part of first term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#' @example PairwiseOneOne12(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3, 3, 2)
PairwiseOneOne12 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B2){
  return((1+(2*kappa+x1+x2)/beta)^{-alpha*rho*B2-1})
}

#' Computes third part of first term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param x2 Positive value corresponding to \code{t2}.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B3 Difference area between \code{t2} and \code{t1} (in this order).
#'
#' @return Second part of first term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#' @example PairwiseOneOne13(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3, 3, 2)
pairwise_11_1_3 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B3){
  return((1+(kappa+x2)/beta)^{-alpha*rho*B3-1})
}

#' Computes first term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param x2 Positive value corresponding to \code{t2}.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B1 Difference area between \code{t1} and \code{t2} (in this order).
#' @param B2 Intersection area between \code{t1} and \code{t2} (in this order).
#' @param B3 Difference area between \code{t2} and \code{t1} (in this order).
#'
#' @return First term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#' @example PairwiseOneOne1(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3, 3, 2)
PairwiseOneOne1 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2, B3){
  temp <- PairwiseOneOne11(t1, x1, t2, x2, alpha, beta, kappa, rho, B1)
  temp <- temp * PairwiseOneOne12(t1, x1, t2, x2, alpha, beta, kappa, rho, B2)
  temp <- temp * pairwise_11_1_3(t1, x1, t2, x2, alpha, beta, kappa, rho, B3)
  return(temp)
}

#' See PairwiseOneOne2 for description.
pairwise_11_2_1 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2){
  return(B1*B2*(1+(2*kappa+x1+x2)/beta)*(1+(kappa+x2)/beta))
}

#' See PairwiseOneOne2 for description.
pairwise_11_2_2 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B3){
  return(B1*B3*(1+(2*kappa+x1+x2)/beta)^2)
}

#' See PairwiseOneOne2 for description.
pairwise_11_2_3 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B2){
  temp <- B2*(B2+1/(alpha*rho))
  temp <- temp*(1+(kappa+x1)/beta)*(1+(kappa+x2)/beta)
  return(temp)
}

#' See PairwiseOneOne2 for description.
pairwise_11_2_4 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B2, B3){
  return(B2*B3*(1+(kappa+x1)/beta)*(1+(2*kappa+x1+x2)/beta))
}

#' Computes second term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param x2 Positive value corresponding to \code{t2}.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param B1 Difference area between \code{t1} and \code{t2} (in this order).
#' @param B2 Intersection area between \code{t1} and \code{t2} (in this order).
#' @param B3 Difference area between \code{t2} and \code{t1} (in this order).
#'
#' @return Second term in latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#' @example PairwiseOneOne2(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, 3, 3, 2)
PairwiseOneOne2 <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2, B3){
  temp <- pairwise_11_2_1(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B2)
  temp <- temp + pairwise_11_2_2(t1, x1, t2, x2, alpha, beta, kappa, rho, B1, B3)
  temp <- temp + pairwise_11_2_3(t1, x1, t2, x2, alpha, beta, kappa, rho, B2)
  temp <- temp + pairwise_11_2_4(t1, x1, t2, x2, alpha, beta, kappa, rho, B2, B3)
  return(temp)
}

#' Computes latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}
#' with exponential trawl function.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param x2 Positive value corresponding to \code{t2}.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param transformation Boolean to use the Marginal Transform (MT) method.
#' @param n_moments Number of moments achieved by transformed GPD marginals, if used.
#'
#' @return Latent trawl pairwise likelihood with \code{(x,x)} where \code{x > 0} and \code{y > 0}.
#' @example PairwiseOneOne2(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, T, 4)
pairwise_11_exp <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation=F, n_moments=4){
  # Marginal Transformation
  if(transformation){
    offset_shape <- n_moments + 1
    offset_scale <- TrfFindOffsetScale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
    inv_x1 <- TrfInverseG(x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    inv_x2 <- TrfInverseG(x2, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    new_x1 <- inv_x1
    new_x2 <- inv_x2
    jacobian1 <- TrfJacobian(z = x1, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    jacobian2 <- TrfJacobian(z = x2, alpha = alpha, beta = beta, kappa = kappa, offset_scale = offset_scale, offset_shape = offset_shape)
    temp <- jacobian1 * jacobian2
    #temp <- 1.0
    #temp <- 1/temp
  }else{
    new_x1 <- x1
    new_x2 <- x2
    temp <- 1.0
  }

  if(temp == 0.0 || is.na(temp) || is.nan(temp)){
    temp <- 1e-16
  }

  trawlA <- ComputeAExp(rho)
  B1 <- ComputeB1Exp(rho, t1, t2)
  B2 <- ComputeBInterExp(rho, t1, t2)
  B3 <- ComputeB3Exp(rho, t1, t2)

  if(transformation){
    temp <- temp * PairwiseOneOne1(t1, new_x1, t2, new_x2, alpha = offset_shape, beta = 1, kappa, rho, B1, B2, B3)
    temp <- temp * PairwiseOneOne2(t1, new_x1, t2, new_x2, alpha = offset_shape, beta = 1, kappa, rho, B1, B2, B3)
  }else{
    temp <- temp * PairwiseOneOne1(t1, new_x1, t2, new_x2, alpha, beta, kappa, rho, B1, B2, B3)
    temp <- temp * PairwiseOneOne2(t1, new_x1, t2, new_x2, alpha, beta, kappa, rho, B1, B2, B3)
  }

  return(temp)
}

# Example
t1 <- 0.0
t2 <- 1.0
x1 <- 1.0
x2 <- 2.0
alpha <- -1.0
beta <- 10
rho <- 1.0
kappa <- 1.0
B1 <- ComputeB1Exp(rho, t1, t2)
B2 <- ComputeBInterExp(rho, t1, t2)
B3 <- ComputeB3Exp(rho, t1, t2)

pairwise_11_exp(t1 = t1, x1 = x1,
                t2 = t2, x2 = x2,
                alpha = alpha, beta = beta,
                rho = rho, kappa = kappa)

answer <- 1/100*(1+0.2)^{B1-1}*(1+0.5)^{B2-1}*(1+3/10)^{B3-1}
answer <- answer * (B1*B2*(1+0.5)*(1+0.3)+B1*B3*(1+0.5)^2
  +B2*(B2-1)*(1+0.2)*(1+0.3)+B2*B3*(1+0.2)*(1+0.5))
answer

#' Computes correct latent trawl SINGLE pairwise likelihood depending on the values of \code{(x1,x2)}
#' with exponential trawl function.
#'
#' @param t1 First timestamp.
#' @param x1 Positive value corresponding to \code{t1}.
#' @param t2 Second timestamp.
#' @param x2 Positive value corresponding to \code{t2}.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param transformation Boolean to use the Marginal Transform (MT) method.
#'
#' @return SINGLE latent trawl pairwise likelihood depending on \code{(x1,x2)}.
#' @example SinglePairPL(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, F)
SinglePairPL <- function(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation=F){
  # TODO check whether t1 should be <= t2 or not
  #print(x1)
  if(x1 < 1e-16){
    if(x2 < 1e-16){
      return(PairwiseZeroZeroExp(t1, t2, alpha = 1, beta = 1, kappa, rho))
    }else{
      return(PairwiseOneZeroExp(t2, x2, t1, alpha, beta, kappa, rho, transformation, n_moments = 0))
    }
  }else{
    if(x2 < 1e-16){
      return(PairwiseOneZeroExp(t1, x1, t2, alpha, beta, kappa, rho, transformation, n_moments = 0))
    }else{
      return(pairwise_11_exp(t1, x1, t2, x2, alpha, beta, kappa, rho, transformation, n_moments = 0))
    }
  }
}

# Example
SinglePairPL(0.1, 2.0, 0.3, 5.0, 2., 3., 30, 0.3, F)
SinglePairPL(0.1, 2.0, 0.3, 1.0, -2., 3., 30, 0.3, T)

#' Computes latent trawl FULL pairwise likelihood depending with exponential trawl function.
#'
#' @param times Vector of timestamps.
#' @param values Vector of target values.
#' @param alpha Shape parameter. Should be positive.
#' @param beta Latent Gamma scale parameter. Should be positive.
#' @param kappa Exceedance probability parameter. Should be positive.
#' @param rho Exponential trawl parameter. Should be positive.
#' @param delta Maximum depth of pairwise likelihood blocks. Should be positive natural integer.
#' @param logscale Boolean to use logscale (log-likelihood). Default \code{T}.
#' @param transformation Boolean to use the Marginal Transform (MT) method.
#'
#' @return FULL latent trawl pairwise likelihood.
#' @example FullPL(t1=1, x1=0.5, t2=4, 0.3, 2, 3, 0.2, F)
FullPL <- function(times, values, alpha, beta, kappa, rho, delta, logscale=T, transformation=F){
  ok_ind <- which(!is.na(values))
  values <- values[ok_ind]
  times <- times[ok_ind]
  k <- length(values)

  temp <- 0.0
  upper <- pmin(1:(k-1)+delta, k)
  lower <- 2:k

  for(i in 1:(k-1)){
    ind <- (lower[i]):(upper[i])
    m <- 0
    for(j in ind){
      warnon <- SinglePairPL(times[i], values[i],
                              times[j], values[j],
                              alpha = alpha,
                              beta = beta,
                              kappa = kappa,
                              rho = rho,
                              transformation=transformation)
      if(!is.na(warnon) & !is.nan(warnon)){
        if(warnon > 1e-12){
          # log the result
          temp <- temp + log(warnon)
        }else{
          if(warnon >= 0.0){
            temp <- temp + log(warnon)
          }
        }
      }
    }
  }

  if(logscale){
    return(temp)
  }else{
    return(exp(temp))
  }
}

#' Wrapper for FullPL using a list for parameters in the exponential trawl case.
#' @param times Vector of timestamps.
#' @param values Vector of target values.
#' @param delta Maximum depth of pairwise likelihood blocks. Should be positive natural integer.
#' @param params List of parameters.
#' @param logscale Boolean to use logscale (log-likelihood). Default \code{T}.
#' @param transformation Boolean to use the Marginal Transform (MT) method.
#'
#' @return Pairwise Likelihood as per FullPL using a list of parameters instead.
#' @example ParamsListFullPL(c(1,2,3,4,5), c(0,2.3,.3,0,5), delta=2,
#' params=list("alpha"=2,"beta"=3,"kappa"=1.5), T, F)
ParamsListFullPL <- function(times, values, delta, params, logscale=T, transformation=F){
  # TODO add general model parameter names
  #print(params)
  return(FullPL(times, values,
                 alpha = (params["alpha"][[1]]),
                 beta = (params["beta"][[1]]),
                 kappa = params["kappa"][[1]],
                 rho = params["rho"][[1]],
                 delta = delta,
                 logscale = T,
                 transformation = transformation))
}

#' Computes univariate latent trawl FULL pairwise likelihood depending with exponential trawl function
#' with the option to fix some parameter values.
#'
#' @param times Vector of timestamps.
#' @param values Vector of target values.
#' @param delta Maximum depth of pairwise likelihood blocks. Should be positive natural integer.
#' @param fixed_names Vector of literal names of parameters to keep fixed.
#' @param fixed_params Vector of numerical values of fixed parameters.
#' @param params List of parameters.
#' @param model_vars_names Vector of all parameters names in the model.
#' @param logscale Boolean to use logscale (log-likelihood). Default \code{T}.
#' @param transformation Boolean to use the Marginal Transform (MT) method.
#'
#' @return FULL latent trawl pairwise likelihood with some (or none) parameters fixed.
#' @example
#' times <- c(1,2,3,4,5)
#' values <- c(2,0,3,4,0)
#' delta <- 2
#' fixed_names <- c("alpha", "beta")
#' params <- c(2.0, 3.4, 0.1, 4.3)
#' model_vars_names <- c("alpha", "beta", "rho", "kappa")
#' UnivariateFullPL(times, values, delta, fixed_names, params, model_vars_names, T, F)
UnivariateFullPL <- function(times, values, delta, fixed_names, fixed_params, params, model_vars_names, logscale=T, transformation=F){
  if(length(fixed_names) > length(model_vars_names)) stop('Too many fixed parameters compared to number of model params.')
  if(length(fixed_params) + length(params) != length(model_vars_names)) stop('Wrong number of params compared to model specs.')
  if(length(fixed_params) != length(fixed_names)) stop('fixed_params and fixed_names should have same length.')

  opti_params <- !(model_vars_names %in% fixed_names)
  opti_params_names <- model_vars_names[opti_params]
  params_all <- rep(0, length(model_vars_names))
  params_all[opti_params] <- params

  if(length(fixed_params) > 0){
    params_all[!opti_params] <- fixed_params
  }

  params_list <- list()
  params_list[fixed_names] <- fixed_params
  #params_list[opti_params_names] <- params
  for(i in 1:length(fixed_names)){
    params_list[fixed_names[i]] <- fixed_params[i]
  }

  params_list[opti_params_names] <- params

  return(ParamsListFullPL(times = times,
                          values = values,
                          delta = delta,
                          params = params_list,
                          logscale = logscale,
                          transformation = transformation))
}

#' Computes Generalised Pareto (log-)likelihood on non-zero exceedances under independence.
#'
#' @param values Vector of target values.
#' @param fixed_names Vector of literal names of parameters to keep fixed.
#' @param fixed_params Vector of numerical values of fixed parameters.
#' @param params List of parameters.
#' @param model_vars_names Vector of all parameters names in the model.
#' @param logscale Boolean to use logscale (log-likelihood). Default \code{T}.
#' @param transformation Boolean to use the Marginal Transform (MT) method.
#' @param n_moments Number of moments the transformed variables should have
#' using the Marginal Transform (MT) method.
#'
#' @return Generalised Pareto (log-)likelihood on non-zero exceedances under independence.
#' @example
#' times <- c(1,2,3,4,5)
#' values <- c(2,0,3,4,0)
#' delta <- 2
#' fixed_names <- c("alpha", "kappa")
#' params <- c(2.0, 3.4, 0.1, 4.3)
#' model_vars_names <- c("alpha", "beta", "rho", "kappa")
#' UnivariateFullPL(times, values, delta, fixed_names, params, model_vars_names, T, F)
MarginalGPDLikelihood <- function(values, fixed_names, fixed_params, params, model_vars_names, logscale=T, transformation=F, n_moments=4){
  if(length(fixed_names) > length(model_vars_names)) stop('Too many fixed parameters compared to number of model params.')
  if(length(fixed_params) + length(params) != length(model_vars_names)) stop('Wrong number of params compared to model specs.')
  if(length(fixed_params) != length(fixed_names)) stop('fixed_params and fixed_names should have same length.')

  opti_params <- !(model_vars_names %in% fixed_names)
  opti_params_names <- model_vars_names[opti_params]
  params_all <- rep(0, length(model_vars_names))
  params_all[opti_params] <- params

  if(length(fixed_params) > 0){
    params_all[!opti_params] <- fixed_params
  }

  if(transformation){
    if(length(params_all) < 3) stop('Marginal GPD with transformation requires 3 parameters: alpha, beta and kappa.')
    print(params_all)
    lik <- vapply(values,
                  function(x){
                    temp_alpha <- params_all[1]
                    temp_beta <- params_all[2]
                    temp_kappa <- params_all[3]
                    offset_shape <- n_moments + 1
                    offset_scale <- TrfFindOffsetScale(alpha = temp_alpha, beta = temp_beta,
                                                          kappa = temp_kappa, offset_shape = offset_shape)
                    inv_x1 <- TrfInverseG(x1, alpha = temp_alpha, beta = temp_beta, kappa = temp_kappa,
                                        offset_scale = offset_scale, offset_shape = offset_shape)
                    jacobian1 <- TrfJacobian(z = x1, alpha = temp_alpha, beta = temp_beta, kappa = temp_kappa,
                                              offset_scale = offset_scale, offset_shape = offset_shape)
                    temp <- jacobian1
                      return(dlgpd(x = inv_x1, alpha = offset_shape, beta = offset_scale+temp_kappa)*jacobian1)
                    },
                  1.0)
  }else{
    lik <- vapply(values, function(x){return(dlgpd(x = x, alpha = params_all[1], beta = (params_all[2]+params_all[3])))}, 1.0)
  }

  if(logscale){
    return(sum(log(lik[lik > 0.0])))
  }else{
    return(prod(lik))
  }
}

#' Simplified marginal (GPD) log-likelihood function under independence in the
#' exponential trawl case.
#'
#' @param values Exceedance values.
#' @param params (at least) 2-d vector for GPD parameters
#'
#' @return Log-likelihood of GPD distribued variables (i.e. non-zero exceedances).
#' @example MarginalSimpleLik(c(2.0, 0.3, 6.15, 0, 0.31), c(2.1, 1.17, 0.52, 4.17)) # for GPD(2.1, 1.17)
MarginalSimpleLik <- function(values, params){
  alpha <- params[1]
  beta_kappa <- params[2]
  lik <- alpha / beta_kappa * (1+sign(alpha) * values / beta_kappa)^{-alpha-1}
  return(sum(log(lik[lik>0])))
}

#' Generalised Pareto likelihood maximisation using L-BFGS-B optimisation routine.
#'
#' @param values Exceedance values.
#' @param initial_guess (at least 2-d) Vector for GPD parameters starting values.
#' @param lower Vector of lower bounds limits for optimisation procedure. Default \code{c(0.1, 0.1)}.
#' @param upper Vector of upper bounds limits for optimisation procedure. Default \code{c(20, 20)}.
#'
#' @return Parameters of Log-likelihood maximisation of GPD distribued variables (i.e. non-zero exceedances).
#' @example MarginalSimpleLik(c(2.0, 0.3, 6.15, 0, 0.31), c(2.1, 1.17, 0.52, 4.17))
GPDFit <- function(values, initial_guess, lower=c(0.1, 0.1), upper=c(20, 20)){
  fn_to_optim <- function(x){return(-MarginalSimpleLik(values = values, params = x))}
  ff <- optim(fn_to_optim, par = initial_guess, method = "L-BFGS-B", lower = lower, upper = upper)
  return(ff$par)
}

#' Method of moments (MoM) on a one or multiple exceedance timeseries for the
#' latent trawl model using marginal GPD properties and probability of exceedance.
#'
#' @param values_array Matrix of exceedance timeseries
#'
#' @return Parameters given by a second-order method of moments as well as standard deviation across
#' timeseries for each individual parameter.
#' @example
#' exceed1 <- c(0.1, 0, 0.2, 0, 0, 0, 0.6, 1.5)
#' exceed2 <- c(0, 0.3, 5.2, 0, 0, 3.0, 0, 2.2)
#' val_array <- cbind(exceed1, exceed2)
#' MoMGPD(val_array)
MoMGPD <- function(values_array){
  # workds under the assumption that alpha > 2
  # values_array contains the time series with first axis as time and second as # of time series
  n_dim <- length(values_array[1,])
  n_values <- length(values_array[,1])
  alphas_mom <- rep(0, n_dim)
  betas_mom <- rep(0, n_dim)
  kappas_mom <- rep(0, n_dim)

  for(index in 1:n_dim){
    var_mom <- var(values_array[,index][values_array[,index]>0])
    mean_mom <- mean(values_array[,index][values_array[,index]>0])
    p_mom <- length(values_array[,index][values_array[,index]>0])/n_values

    alphas_mom[index] <- 2*var_mom/(var_mom-mean_mom^2)
    betas_mom[index] <- mean_mom*(alphas_mom[index]-1)

    kappas_mom[index] <- betas_mom[index] * (1.0 - p_mom^{1/alphas_mom[index]})
    betas_mom[index] <- betas_mom[index] - kappas_mom[index]
  }

  return(list(alpha=alphas_mom, beta=betas_mom, kappa=kappas_mom,
              mean_alpha=mean(alphas_mom), mean_beta=mean(betas_mom), mean_kappa=mean(kappas_mom),
              sd_alpha=sd(alphas_mom), sd_beta=sd(betas_mom), sd_kappa=sd(kappas_mom)))
}

