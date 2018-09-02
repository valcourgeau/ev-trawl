#' Sample Gamma with small shape parameter \code{alpha*dx*dy}.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param dx x-axis multiplier.
#' @param dy y-axis multiplier.
#' @param n Sample size.
#'
#' @return n Gamma(alpha*dx*dy, beta) samples.
#' @example GammaBox(2.0, 3.0, 0.5, 0.5, 10)
#' @importFrom stats rgamma
GammaBox <- function(alpha, beta, dx, dy, n){
  requireNamespace("stats", quietly = T)
  alpha <- abs(alpha)
  beta <- abs(beta)
  return(stats::rgamma(shape=alpha*dx*dy, rate=beta, n=n))
}

#' Wrapper for GammaBox using a square to multiply the shape parameter.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param dt Both x-axis and y-axis multiplier.
#' @param n Sample size.
#'
#' @return n Gamma(alpha*dt*dt, beta) samples.
#' @example GammaSqBox(2.0, 3.0, 0.5, 10)
GammaSqBox <- function(alpha, beta, dt, n){
  return(GammaBox(alpha = alpha, beta = beta, dx = dt, dy = dt, n = n))
}

#' Returns an exponential trawl function with base time t.
#'
#' @param t trawl function peak time.
#' @param rho exponential trawl parameter. Should be positive.
#' @param max_value Miximal value of the trawl function (if known). Default 1.
#' @param min_value Minimal value accepted for the trawl function. Default is
#'   1e-2.
#'
#' @return (Vectorised) Exponential trawl function with peak time t and
#'   parameter rho. If this function is evaluated using NA, it yields a list of
#'   key components (rho, max time difference, total area of trawl set A).
#' @example TrawlExp(t=1, rho=1) f(0) # returns e^{-1}
TrawlExp <- function(t, rho, max_value=1, min_value=1e-2){
  if(rho <= 0) stop('rho should be positive.')
  time_eta <- -log(min_value)/rho
  return(function(s){
    if(is.na(s[1])){
      return(list(rho=rho,
                  time_eta=time_eta,
                  max_value=max_value,
                  trawl_time=t,
                  A=ComputeAExp(rho)))
    }

    s_to_use <- which(s <= t & s >= t-time_eta)
    results <- rep(0, length(s))
    results[s_to_use] <- exp(-rho * (t-s[s_to_use]))
    return(results)
  })
}

#' Returns a primitive of exponential trawl function with base time t. Such
#' primitive has its zero at \code{zero_at}.
#'
#' @param t trawl function peak time.
#' @param rho exponential trawl parameter. Should be positive.
#' @param zero_at Value at which primitive is zero. Default is \code{-Inf}.
#'
#' @return (Vectorised) Primitive of Exponential trawl function with peak time
#'   t and parameter rho. If this function is evaluated using NA, it yields a
#'   list of key components such as the trawl peak time t.
#' @example TrawlExpPrimitive(t=1, rho=2) F(0) # returns e^{-2}/2
TrawlExpPrimitive <- function(t, rho, zero_at=-Inf){
  # primitive of exponential trawl function which is zero at zero_at
  if(rho <= 0) stop('rho should be positive.')
  return(function(s){
    if(is.na(s[1])){
      return(list(trawl_time=t))
    }

    s_to_use <- which(s < t)
    results <- rep(1/rho - exp(-rho * (t-zero_at)) / rho, length(s))
    results[s_to_use] <- exp(-rho * (t-s)) / rho - exp(-rho * (t-zero_at)) / rho
    return(results)
  })
}

#' Creates a list of trawl functions  given a type of trawl and a vector of
#' times as peak times. Has the option to get primitives. Note that only
#' exponential trawl are implemented so far.
#'
#' @param times Vector of times to evaluate the trawl function (or primitive)
#'   at.
#' @param params List of trawl parameters.
#' @param type Trawl type (so far, only "exp" is available).
#' @param prim Boolean to use primitive or not. Default is False (F).
#'
#' @return Collection of trawl functions set on \code{times} given the type of
#'   trawl (\code{type}).
#' @example CollectionTrawl(times = 1:5, list("rho"=0.4), type="exp")
CollectionTrawl <- function(times, params, type, prim=F){
  if(!is.list(params)) stop('params should be a list.')
  # TODO Add more than exp
  if(type=="exp"){
    if(! "rho" %in% names(params)) stop('rho should in list of parameters (params).')
    params_rho <- params$rho
    results <-
      for(time_index in 1:length(times)){
        if(prim){
          return(lapply(times, function(t){TrawlExpPrimitive(t, params_rho)}))
        }else{
          return(lapply(times, function(t){TrawlExp(t, params_rho)}))
        }
      }
  }
  stop('Fatal error: no trawl functions list created')
}

#' Computes the area of a slice for the Slice Partition method of trawl
#' functions.
#'
#' @param i main index.
#' @param j secondary index.
#' @param times vector of discret times at which the trawl function is
#'   partionned.
#' @param trawl.f.prim Trawl function primitive.
#'
#' @return Area of slice \code{S(i,j)} in Slice Partition method for trawl
#'   functions.
#' @examples SliceArea(1, 2, times=c(0, 1, 2, 3), TrawlExpPrimitive(t=2,
#'   rho=0.2))
SliceArea <- function(i, j, times, trawl.f.prim){
  prim_info <- trawl.f.prim(NA)
  origin_time <- prim_info$trawl_time
  # TODO make sure times are sorted before using this function
  # times <- sort(times)
  if(i != 1 & j != length(times)){
    temp <- trawl.f.prim(times[i] - times[j] + origin_time)
    temp <- temp - trawl.f.prim(times[i] - times[j+1] + origin_time)
    temp <- temp - trawl.f.prim(times[i-1] - times[j] + origin_time)
    temp <- temp + trawl.f.prim(times[i-1] - times[j+1] + origin_time)
  }else{
    temp <- trawl.f.prim(times[i] - times[j] + origin_time)
    if(j != length(times)){
      temp <- temp - trawl.f.prim(times[i] - times[j+1] + origin_time)
    }else{
      if(i == length(times)){
        temp <- trawl.f.prim(origin_time) - trawl.f.prim(origin_time-times[i]+times[i-1])
      }
    }
  }

  return(temp)
}

#' Performs trawl slices reconstuction to get gamma samples using the
#' independent measure scaterring.
#'
#' @param alpha (Gamma) Shape parameter.
#' @param beta (Gamma) Scale parameter.
#' @param times Vectors of discret times.
#' @param marg.dist Name of infinitely divisible distribution. Currenlty
#'   implemented: gamma, gaussian, generalised hyperbolic (ghyp), generalised
#'   inverse gaussian (gig).
#' @param n Number of simulations (so far, only \code{n=1} is implemented).
#' @param trawl.fs collection of trawl functions indexed on \code{times}.
#' @param trawl.fs.prim collection of trawl functions primitives indexed on
#'   \code{times}.
#' @param deep_cols Depth of reconstruction (columns). Default is 30.
#'
#' @return Samples using trawl slice reconstruction disributed using
#'   \code{marginal}.
#' @examples fcts <- lapply(c(1,2,3,4), function(t) TrawlExp(t, rho=0.2)) prims
#'   <- lapply(c(1,2,3,4), function(t) TrawlExpPrimitive(t, rho=0.2))
#'   TrawlSliceReconstruct(alpha=3, beta=2, times=c(0, 1, 2, 3), marg.dist =
#'   "gamma", n=1, fcts, prims)
#'   @importFrom ghyp ghyp
#'   @importFrom stats rgamma rnorm
TrawlSliceReconstruct <- function(alpha, beta, times, marg.dist, n, trawl.fs, trawl.fs.prim, deep_cols=30, ghyp.object=NA){
  # TODO Add GIG compatibility
  # TODO sort the trawl.fs and trawl.fs.prim as the times
  requireNamespace("ghyp", quietly = TRUE)
  requireNamespace("stats", quietly = TRUE)

  if(n > 1) stop("Case n>1 not yet implemented.")
  if(!is.list(trawl.fs.prim)) stop('Wrong type: trawl function primitives should be a list.')
  if(!marg.dist %in% c("gamma", "normal", "gaussian", "gig", "ghyp")){
    stop(paste('marg.distr', marg.dist, 'not yet implemented.'))
  }else if(marg.dist == "ghyp" & class(ghyp::ghyp())[1] != "ghyp"){
    stop('ghyp.object should be an instance of ghyp::ghyp')
  }

  n_times <- length(times)

  A <- trawl.fs[[1]](NA)$A # TODO A special for each timestep
  slice_mat <- matrix(0, nrow = n_times, ncol = deep_cols)
  gamma_sim <- matrix(0, nrow = n_times, ncol = deep_cols)

  # Creating the matrix of gamma realisations
  for(main_index in 1:n_times){
    for(second_index in 1:deep_cols){
      slice_mat[main_index, second_index] <- SliceArea(main_index, min(second_index + main_index - 1, n_times),
                                                       times, trawl.fs.prim[[main_index]]) # TODO fix for last row

      # it suffices to implement new marginals here
      gamma_sim[main_index, second_index] <-
        switch(marg.dist,
             "gamma" = stats::rgamma(n = 1, alpha * slice_mat[main_index, second_index] / A,
                                beta),
             "gaussian" = stats::rnorm(n = 1, mean = alpha * slice_mat[main_index, second_index] / A,
                                  sd =  beta * sqrt(slice_mat[main_index, second_index] / A)),
             "normal" = stats::rnorm(n = 1, mean = alpha * slice_mat[main_index, second_index] / A,
                              sd = beta * sqrt(slice_mat[main_index, second_index] / A)),
             "gig" = ghyp::rgig(n = 1, lambda = -0.5, chi = alpha*slice_mat[main_index, second_index] / A,
                                psi = beta*slice_mat[main_index, second_index] / A), # TODO
             "gh" = ghyp::rghyp(n=1, object = ghyp.object)
             # TODO find what kind of new variables we need for ID property
             )

    }
  }

  # Using independent scattering of Levy basis to add time dependence to trawls
  results <- matrix(0, nrow = length(times), ncol = 1)
  upper.anti.tri<-function(m){col(m)+row(m) < dim(m)[1]+1}
  anti_diag_mat <- matrix(1, deep_cols, deep_cols)
  anti_diag_mat[upper.anti.tri(anti_diag_mat)] <- 0
  results <- vapply(1:(n_times - 1*deep_cols),
                    function(i){
                      return(sum(gamma_sim[i:(i - 1 + deep_cols),] * anti_diag_mat))
                    },
                    1.0)

  return(results)
}

#' Simulation of trawl process path using Slice partition.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param rho Trawl parameters. Must be positive if Exponential trawls are used.
#' @param times Vectors of discret times.
#' @param marg.dist Name of infinitely divisible distribution. Currenlty
#'   implemented: gamma, gaussian, generalised hyperbolic, generalised inverse
#'   gaussian.
#' @param trawl.function Type of trawl function that should be used. Default NA.
#' @param trawl.fs collection of trawl functions indexed on \code{times}.
#'   Default NA. Default NA if no \code{trawl.function} is indicated and should
#'   contain as many as in \code{times}.
#' @param trawl.fs.prim collection of trawl functions primitives indexed on
#'   \code{times}. Default NA if no \code{trawl.function} is indicated and
#'   should contain as many as in \code{times}.
#' @param n Number of simulations (so far, only \code{n=1} is implemented).
#' @param kappa Additive constant to scale parameter \code{beta}.
#' @param transformation Boolean to apply marginal transform method. Default is
#'   False (F).
#' @param offset_shape Transformed-marginal Shape parameter.
#' @param offset_scale Transformed-marginal Scale parameter.
#' @param deep_cols Depth of reconstruction (columns). Default is 30.
#'
#' @return Simulated path (size the same as times) of trawl process.
#' @examples
#' alpha <- 5
#' beta <- 3
#' times <- 1:15
#' trawl.function <- "exp"
#' margi <- rgamma
#' kappa <- 0
#' rtrawl(alpha = alpha, beta = beta, kappa = kappa, marginals = margi,
#'  trawl.function = trawl.function)
#'
#' @export
rtrawl <- function(alpha, beta, times, marg.dist, trawl.function=NA, trawl.fs=NA, trawl.fs.prim=NA, n, rho=NA,
                   kappa = 0, transformation=F, offset_shape=NULL, offset_scale=NULL, deep_cols=30){
  if(!is.na(trawl.function)){
    if(trawl.function %in% c("exp")){
      if(is.na(rho)) stop('If trawl.function is not NA, need trawl parameters rho.')
      if(rho <= 0) stop('rho should be positive.')

      trawl.fs <- CollectionTrawl(times = times,
                                  params = list("rho"=rho), type = trawl.function)
      trawl.fs.prim <- CollectionTrawl(times = times,
                                       params = list("rho"=rho), type = trawl.function,
                                       prim = T)
    }else{
      stop(paste("trawl.function", trawl.function, "not yet implemented."))
    }

  }else{
    if(length(trawl.fs) != length(times)){
      stop('Wrong number of trawl functions compared to timestamps.')
    }
    if(length(trawl.fs.prim) != length(times)){
      stop('Wrong number of trawl primitives compared to timestamps.')
    }
  }

  if(transformation & (is.null(offset_scale) | is.null(offset_shape))){
    stop('When using marginal trf, indicate shape and scale offsets.')
  }

  if(!transformation){
    results <- TrawlSliceReconstruct(alpha = alpha,
                                      beta = beta+kappa,
                                      marg.dist = marg.dist,
                                      times = times,
                                      trawl.fs = trawl.fs,
                                      trawl.fs.prim = trawl.fs.prim,
                                      n = n, deep_cols)
  }else{
    results <- TrawlSliceReconstruct(alpha = offset_shape,
                                      beta = offset_scale,
                                      marg.dist = marg.dist,
                                      times = times,
                                      trawl.fs = trawl.fs,
                                      trawl.fs.prim = trawl.fs.prim,
                                      n = n, deep_cols)
  }

  return(results)
}

#' Simulation of extreme value path using latent trawl process. Transformed
#' marginals have scale parameter \code{1+kappa}.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param kappa Additive constant to scale parameter \code{beta}.
#' @param rho Trawl parameters. Must be positive if Exponential trawls are used.
#' @param times Vectors of discret times.
#' @param marg.dist Name of infinitely divisible distribution for latent trawls.
#'   Currenlty implemented: gamma, gaussian, generalised hyperbolic, generalised
#'   inverse gaussian.
#' @param n Number of simulations (so far, only \code{n=1} is implemented).
#' @param transformation Boolean to apply marginal transform method. Default is
#'   False (F).
#' @param trawl.function Type of trawl function that should be used. Default NA.
#' @param trawl.fs collection of trawl functions indexed on \code{times}.
#' @param trawl.fs.prim collection of trawl functions primitives indexed on
#'   \code{times}.
#' @param n_moments Number of finite moments for transformed marginals.
#' @param deep_cols Depth of reconstruction (columns). Default is 30.
#'
#' @return Simulated path (size the same as times) of latent-trawl extreme
#'   value process.
#' @examples
#' alpha <- 3
#' beta <- 2
#' kappa <- 0.95
#' rho <- 0.2
#' n.timestamps <- 200
#' times <- 1:n.timestamps
#'
#' marg.dist <- "gamma"
#' n <- 1
#' transformation <- F
#' trawl.function <- "exp"
#'
#' rlexceed(alpha = alpha, beta = beta, kappa = kappa, rho = rho, times = times,
#'          marg.dist = marg.dist, n = n, transformation = transformation,
#'          trawl.function= trawl.function)
#'
#' @importFrom stats runif rexp
#' @export
rlexceed <- function(alpha, beta, kappa, rho=NA, times, marg.dist, n, transformation,
                     trawl.function=NA, trawl.fs=NA, trawl.fs.prim=NA, n_moments = 4, deep_cols=30){
  requireNamespace("stats", quietly = T)

  # TODO n > 1 is not implemented yet
  offset_shape <- n_moments+1
  offset_scale <- TrfFindOffsetScale(alpha = alpha,
                                        beta = beta,
                                        kappa = kappa,
                                        offset_shape = n_moments+1)
  #print(offset_scale)
  # Generate latent process
  gen_trawl <- rtrawl(alpha = alpha,
                       beta = beta,
                       kappa = 0.0,
                       rho = rho,
                       marg.dist = marg.dist,
                       times = times,
                       trawl.function = trawl.function,
                       trawl.fs = trawl.fs,
                       trawl.fs.prim = trawl.fs.prim,
                       n = n,
                       transformation = transformation,
                       offset_shape = offset_shape,
                       offset_scale = 1)

  #return(gen_trawl)
  # Uniform threshold
  unif_samples <- stats::runif(n=length(times)-deep_cols)
  if(n == 1){
    gen_exceedances <- rep(0, length(times)-deep_cols)
  }else{
    gen_exceedances <- matrix(0, nrow = length(times)-deep_cols, ncol = n)
  }

  #print(gen_trawl)
  prob_zero <- 1.0-exp(-kappa * gen_trawl)
  which_zero <- which(prob_zero >= unif_samples)
  if(transformation){
    gen_exceedances[-which_zero] <-  vapply(stats::rexp(n = length(gen_trawl[-which_zero]), rate = gen_trawl[-which_zero]),
                                            FUN.VALUE = 1.0,
                                            FUN = function(x){return(TrfG(x, alpha = alpha,
                                                                           beta = beta,
                                                                           kappa = kappa,
                                                                           offset_scale = offset_scale,
                                                                           offset_shape = offset_shape))})
  }else{
    gen_exceedances[-which_zero] <-  stats::rexp(n = length(gen_trawl[-which_zero]), rate = gen_trawl[-which_zero])
  }

  return(gen_exceedances)
}
