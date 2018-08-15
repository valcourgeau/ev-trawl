#' Centered differentiation calculation.
#'
#' @param value.p function value at \code{x+epsilon}.
#' @param value.m function value at \code{x-epsilon}.
#' @param epsilon offset value for \code{x}.
#'
#' @return \code{(value.p - value.m) / (2*epsilon)}.
#' @examples
#' DiffVal(4, 2, 1)
#' DiffVal(2.5, 2, 2)
DiffVal <- function(value.p, value.m, epsilon){
  # n

  return((value.p-value.m)/(2*epsilon))
}

#' Second-order centered differentiation calculation.
#'
#' @param value.p function value at \code{x+epsilon}.
#' @param value.c function value at \code{x}.
#' @param value.m function value at \code{x-epsilon}.
#' @param epsilon offset value for \code{x}.
#'
#' @return \code{(value.p -2*value.c + value.m) / (2*epsilon)}.
#' @examples
#' SecondOrderDiffVal(4, 3, 2, 1e-6)
#' SecondOrderDiffVal(2.5, 5, 2, 1e-6)
SecondOrderDiffVal <- function(value.p, value.c, value.m, epsilon){
  temp <- value.p + value.m
  temp <- temp - 2.0*value.c
  temp <- temp / epsilon
  return(temp / epsilon)
}

#' Second-order bivariate differentiation calculation.
#'
#' @param value.c function value at \code{(x,y)}.
#' @param value.x.p function value at \code{(x+epsilonx,y)}.
#' @param value.x.m function value at \code{(x-epsilon,y)}.
#' @param value.x.p function value at \code{(x,y+epsilonx)}.
#' @param value.x.m function value at \code{(x,y-epsilon)}.
#' @param epsilon offset value for both \code{x} and \code{y}.
#'
#' @return \code{(value.p -2*value.c + value.m) / (2*epsilon)}.
#' @examples
#' SecondOrderMixedDiffVal(4, 3, 2, 3.5, 3, 2, 2.12, 1e-6)
#' SecondOrderMixedDiffVal(4.52, 3.14, 2, 3.5, 3.175, 2.527, 2.51, 1e-1)
SecondOrderMixedDiffVal <- function(value.c,
                             value.x.p, value.x.m,
                             value.y.p, value.y.m,
                             value.xy.p, value.xy.m,
                             epsilon){
  temp <- + value.xy.p + value.xy.m  + 2.0*value.c
  temp <- temp - value.x.m - value.y.m - value.x.p - value.y.p
  temp <- temp / epsilon
  return(temp/(2.0*epsilon))
}

#' Computes gradident of a function given parameters and offset value epsilon.
#'
#' @param f R function taking only a vector of parameter as input.
#' @param params Vector of parameters.
#' @param epsilon Offset value for all components.
#'
#' @return Gradient of function \code{f} at parameters \code{params} with offset \code{epsilon}.
#'
#' @examples
#' GradF(function(params){prod(params^2)}, params = c(0.5, 7), epsilon = 1e-6)
#' GradF(function(params){sum(params)}, params = c(2.5, 7), epsilon = 1e-6)
GradF <- function(f,
                   params,
                   epsilon=1e-6){
  params.fixed <- params
  d <- length(params)
  answer <- rep(0, d)

  for(index_par in 1:d){
    params.fixed[index_par] <- params[index_par] + epsilon
    value.p <- f(params.fixed)

    params.fixed[index_par] <- params[index_par] - epsilon
    value.m <- f(params.fixed)

    answer[index_par] <- DiffVal(value.p = value.p,
                                  value.m = value.m,
                                  epsilon = epsilon)
    params.fixed <- params
  }

  return(answer)
}

#' Evalutae function at given parameters.
#'
#' @param f R function taking only a vector of parameters as input.
#' @param params Vector of parameters.
#' @param epsilon Offset value to add to all component of \code{params}.
#'
#' @return Evaluate of function \code{f} at parameters \code{params} with offset \code{epsilon}.
#'
#' @examples
#' EvaluateF(function(params){prod(params^2)}, params = c(0.5, 7), epsilon = 0)
#' EvaluateF(function(params){sum(params)}, params = c(2.5, 7), epsilon = 1)
EvaluateF <- function(f,
                       params,
                       epsilon){
  d <- length(params)
  eval.f <- matrix(0, d, d)
  params.fixed <- params

  for(main in 1:d){
    params.fixed[main] <- params[main] + epsilon
    for(second in main:d){
      if(second > main){
        params.fixed[second] <- params[second] + epsilon
        eval.f[main, second] <- f(params.fixed)
        params.fixed[second] <- params[second]
      }else{
        eval.f[main, second] <- f(params.fixed)
      }

    }

    params.fixed <- params
  }

  return(eval.f)
}

#' Computes Hessian matrix of a $C^2$ function given parameters and offset value epsilon
#'
#' @param f R function taking only a vector of parameter as input.
#' @param params Vector of parameters.
#' @param epsilon Offset value for all components.
#'
#' @return Hessian matrix of function \code{f} at parameters \code{params} with offset \code{epsilon}.
#'
#' @examples
#' HessianF(function(params){prod(params^2)}, params = c(0.5, 7), epsilon = 1e-6)
#' HessianF(function(params){sum(params^3)}, params = c(2.5, 7), epsilon = 1e-6)
HessianF <- function(f,
                      params,
                      epsilon=1e-6){
  if(epsilon < 0.0) stop("Epsilon must be positive.")

  d <- length(params)
  eval.f.p <- EvaluateF(f, params, epsilon = epsilon)
  eval.f.m <- EvaluateF(f, params, epsilon = -epsilon)
  hess.f <- matrix(0.0, d, d)
  f.value <- f(params)

  for(main in 1:d){
    for(second in main:d){
      if(main == second){
        hess.f[main, second] <- SecondOrderDiffVal(value.p = eval.f.p[main, second],
                                             value.m = eval.f.m[main, main],
                                             value.c = f.value,
                                             epsilon = epsilon)
      }else{
        hess.f[main, second] <- SecondOrderMixedDiffVal(value.x.p = eval.f.p[main, main],
                                                 value.x.m = eval.f.m[main, main],
                                                 value.y.p = eval.f.p[second, second],
                                                 value.y.m = eval.f.m[second, second],
                                                 value.xy.p = eval.f.p[main, second],
                                                 value.xy.m = eval.f.m[main, second],
                                                 value.c = f.value,
                                                 epsilon = epsilon)
      }
    }
  }

  hess.f <- hess.f + t(hess.f)
  diag(hess.f) <- 0.5*diag(hess.f)
  return(hess.f)
}
