setwd("~/GitHub/ev.trawl/R/")
source('pairwise_latent_trawl.R')

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
GammaBox <- function(alpha, beta, dx, dy, n){
  alpha <- abs(alpha)
  beta <- abs(beta)
  return(rgamma(shape=alpha*dx*dy, rate=beta, n=n))
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
#' @returns (Vectorised) Exponential trawl function with peak time t and
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
#' @returns (Vectorised) Primitive of Exponential trawl function with peak time
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
#' @returns Collection of trawl functions set on \code{times} given the type of
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
#' @param trawl_f_prim Trawl function primitive.
#'
#' @returns Area of slice \code{S(i,j){} in Slice Partition method for trawl
#'   functions.
#' @example SliceArea(1, 2, times=c(0, 1, 2, 3), TrawlExpPrimitive(t=2,
#'   rho=0.2))
SliceArea <- function(i, j, times, trawl_f_prim){
  prim_info <- trawl_f_prim(NA)
  origin_time <- prim_info$trawl_time
  # TODO make sure times are sorted before using this function
  # times <- sort(times)
  if(i != 1 & j != length(times)){
    temp <- trawl_f_prim(times[i] - times[j] + origin_time)
    temp <- temp - trawl_f_prim(times[i] - times[j+1] + origin_time)
    temp <- temp - trawl_f_prim(times[i-1] - times[j] + origin_time)
    temp <- temp + trawl_f_prim(times[i-1] - times[j+1] + origin_time)
  }else{
    temp <- trawl_f_prim(times[i] - times[j] + origin_time)
    if(j != length(times)){
      temp <- temp - trawl_f_prim(times[i] - times[j+1] + origin_time)
    }else{
      if(i == length(times)){
        temp <- trawl_f_prim(origin_time) - trawl_f_prim(origin_time-times[i]+times[i-1])
      }
    }
  }

  return(temp)
}
#
# trawl_slice_sets_not_optim <- function(alpha, beta, times, n, trawl_fs, trawl_fs_prim){
#   # TODO sort the trawl_fs and trawl_fs_prim as the times
#   # TODO current_trawl and going further back? instead of 8
#
#   if(!is.list(trawl_fs)) stop('Wrong type: trawl set should be a list.')
#
#   times <- sort(times)
#   A <- trawl_fs[[1]](NA)$A
#   slice_mat <- matrix(0, nrow = length(times), ncol = length(times))
#   gamma_sim <- matrix(0, length(times) * length(times), ncol= n)
#
#   # Creating the matrix of gamma realisations
#   for(main_index in 1:length(times)){
#     for(second_index in main_index:length(times)){
#       slice_mat[main_index, second_index] <- SliceArea(main_index, second_index, times, trawl_fs_prim[[main_index]])
#       gamma_sim[(main_index-1) * length(times) + second_index,] <- rgamma(shape = alpha * slice_mat[main_index, second_index] / A,
#                                                                           rate = beta,
#                                                                           n = n)
#     }
#   }
#
#   # Going back in time to use the dep structure
#   n_trawl_forward <- trawl_fs[[1]](NA)$time_eta
#   results <- matrix(0, nrow = length(times), ncol = n)
#   for(current_trawl in 1:length(times)){
#     for(back_trawl in max(1, current_trawl-8+1):current_trawl){
#       for(slice_piece in (current_trawl):min(floor(current_trawl+n_trawl_forward+1), length(times))){
#         if(back_trawl != current_trawl | slice_piece > current_trawl){
#           results[current_trawl,] <- results[current_trawl,] + gamma_sim[(back_trawl-1) * length(times) + slice_piece,]
#         }
#       }
#     }
#   }
#
#   # Using independent scattering of Levy basis to add time dependence to trawls
#   for(main_index in 1:(length(times))){
#     results[main_index,] <- results[main_index,] + gamma_sim[(main_index-1) * length(times)+main_index,]
#   }
#
#   #results[length(times),] <- results[length(times), ] + gamma_sim[(length(times)-1) * length(times) + length(times),]
#
#   return(results)
# }
#
# trawl_slice_sets_not_optim_2 <- function(alpha, beta, times, n, trawl_fs, trawl_fs_prim){
#   # TODO sort the trawl_fs and trawl_fs_prim as the times
#   # TODO current_trawl and going further back? instead of 8
#
#   if(!is.list(trawl_fs)) stop('Wrong type: trawl set should be a list.')
#
#   deep_cols <- 30 # TODO render this customisable
#
#   #times <- sort(times)
#   A <- trawl_fs[[1]](NA)$A # TODO A special for each timestep
#   slice_mat <- matrix(0, nrow = length(times), ncol = deep_cols)
#   gamma_sim <- matrix(0, length(times) * deep_cols, ncol= n)
#
#   # Creating the matrix of gamma realisations
#   for(main_index in 1:length(times)){
#     for(second_index in 1:deep_cols){
#       slice_mat[main_index, second_index] <- SliceArea(main_index, min(second_index + main_index - 1, length(times)), times, trawl_fs_prim[[main_index]])
#       #if(slice_mat[main_index, second_index] > 1e-3){
#       #print(alpha * slice_mat[main_index, second_index] / A)
#       gamma_sim[(main_index-1) * deep_cols + second_index,] <- rgamma(shape = alpha * slice_mat[main_index, second_index] / A,
#                                                                       rate = beta,
#                                                                       n = n)
#       #}
#     }
#   }
#
#   # Going back in time to use the dep structure
#   # n_trawl_forward <- trawl_fs[[1]](NA)$time_eta
#   n_trawl_forward <- deep_cols # TODO remove this
#
#   # Using independent scattering of Levy basis to add time dependence to trawls
#   results <- matrix(0, nrow = length(times), ncol = n)
#   for(current_trawl in 1:length(times)){
#     for(back_trawl in max(1, current_trawl-15+1):current_trawl){
#       for(slice_piece in (current_trawl-back_trawl+1):min(deep_cols, length(times)-back_trawl+1)){
#         if(back_trawl != current_trawl | slice_piece > 1){
#           #if(sum(abs(gamma_sim[(back_trawl-1) * deep_cols + slice_piece,])) > 1e-16){
#           results[current_trawl,] <- results[current_trawl,] + gamma_sim[(back_trawl-1) * deep_cols + slice_piece,]
#           # }
#         }
#       }
#     }
#   }
#
#   # Adding the part that is unique to each column: the top
#   for(main_index in 1:(length(times))){
#     results[main_index,] <- results[main_index,] + gamma_sim[(main_index-1) * deep_cols + 1,]
#   }
#
#   #results[length(times),] <- results[length(times), ] + gamma_sim[(length(times)-1) * length(times) + length(times),]
#
#   return(results)
# }

#' Performs trawl slices reconstuction to get gamma samples using the
#' independent measure scaterring.
#'
#' @param alpha (Gamma) Shape parameter.
#' @param beta (Gamma) Scale parameter.
#' @param times Vectors of discret times.
#' @param marginal Function to generate marginals.
#' @param n Number of simulations (so far, only \code{n=1} is implemented).
#' @param trawl_fs collection of trawl functions indexed on \code{times}.
#' @param trawl_fs_prim collection of trawl functions primitives indexed on
#'   \code{times}.
#' @param deep_cols Depth of reconstruction (columns). Default is 30.
#'
#' @return Samples using trawl slice reconstruction disributed using \code{marginal}.
#' @example fcts <- lapply(c(1,2,3,4), function(t) TrawlExp(t, rho=0.2))
#' prims <- lapply(c(1,2,3,4), function(t) TrawlExpPrimitive(t, rho=0.2))
#' TrawlSliceReconstruct(alpha=3, beta=2, times=c(0, 1, 2, 3), n=1, fcts, prims)
TrawlSliceReconstruct <- function(alpha, beta, times, marginal, n, trawl_fs, trawl_fs_prim, deep_cols=30){
  # TODO Function for other marginals
  # TODO REMOVE trawl_fs dependency
  # TODO check size times vs trawl_fs_prim
  # TODO sort the trawl_fs and trawl_fs_prim as the times
  # TODO current_trawl and going further back? instead of 8
  if(n > 1) stop("Case n>1 not yet implemented.")
  if(!is.list(trawl_fs_prim)) stop('Wrong type: trawl set should be a list.')

  #deep_cols <- 30 # TODO render this customisable
  n_times <- length(times)

  A <- trawl_fs[[1]](NA)$A # TODO A special for each timestep
  slice_mat <- matrix(0, nrow = n_times, ncol = deep_cols)
  gamma_sim <- matrix(0, nrow = n_times, ncol = deep_cols)

  # Creating the matrix of gamma realisations
  for(main_index in 1:n_times){
    for(second_index in 1:deep_cols){
      slice_mat[main_index, second_index] <- SliceArea(main_index, min(second_index + main_index - 1, n_times), times, trawl_fs_prim[[main_index]]) # TODO fix for last row
      gamma_sim[main_index, second_index] <- marginal(shape = alpha * slice_mat[main_index, second_index] / A,
                                                    rate = beta,
                                                    n = 1)
    }
  }

  # Using independent scattering of Levy basis to add time dependence to trawls
  results <- matrix(0, nrow = length(times), ncol = 1)
  upper.anti.tri<-function(m) col(m)+row(m) < dim(m)[1]+1
  anti_diag_mat <- matrix(1, deep_cols, deep_cols)
  anti_diag_mat[upper.anti.tri(anti_diag_mat)] <- 0
  results <- vapply(1:(n_times - 1*deep_cols),
                    function(i){
                      return(sum(gamma_sim[i:(i - 1 + deep_cols),] * anti_diag_mat))
                    },
                    1.0)

  return(results)
}

# ftrawl_slice_sets <- function(alpha, beta, times, trawl_f, trawl_f_prim, n){
#   times <- sort(times)
#   trawl_info <- trawl_f(NA)
#   A <- trawl_info$A
#   results <- matrix(0, nrow = length(times), ncol = n)
#
#   for(i in 1:length(times)){
#     max_index <- suppressWarnings(min(length(times), min(which(times - times[i] > trawl_info$time_eta))))
#     for(j in (i+1):length(times)){
#       if(j <= length(times)){
#         results[i,] <- results[i,] + rgamma(shape = alpha * SliceArea(i, j, times, trawl_f_prim) / A,
#                                             rate = beta,
#                                             n = n)
#       }
#     }
#     if(i < length(times)){
#       results[i+1,] <- results[i+1,] + results[i,]
#       results[i,] <- results[i,] + rgamma(shape = alpha * SliceArea(i, i, times, trawl_f_prim) / A,
#                                           rate = beta,
#                                           n = n)
#     }
#   }
#
#   results[length(times),] <- results[length(times), ] + rgamma(shape = alpha * SliceArea(length(times), length(times), times, trawl_f_prim) / A,
#                                                                rate = beta,
#                                                                n = n)
#   return(results)
# }

#' Simulation of trawl process path using Slice partition.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param times Vectors of discret times.
#' @param marginal Function to generate marginals.
#' @param trawl.function Type of trawl function that should be used. Default NA.
#' @param trawl_fs collection of trawl functions indexed on \code{times}.
#'   Default NA. Default NA if no \code{trawl.function} is indicated and should
#'   contain as many as in \code{times}.
#' @param trawl_fs_prim collection of trawl functions primitives indexed on
#'   \code{times}. Default NA if no \code{trawl.function} is indicated and
#'   should contain as many as in \code{times}.
#' @param n Number of simulations (so far, only \code{n=1} is implemented).
#' @param kappa Additive constant to scale parameter \code{beta}.
#' @param transformation Boolean to apply marginal transform method. Default is
#'   False (F).
#' @param offset_shape Transformed-marginal Shape parameter.
#' @param offset_scale Transformed-marginal Scale parameter.
#'
#' @returns Simulated path (size the same as times) of trawl process.
#' @example
#' alpha <- 5
#' beta <- 3
#' times <- 1:15
#' trawl.function <- "exp"
#' margi <- rgamma
#' kappa <- 0
#' rtrawl(alpha = alpha, beta = beta, kappa = kappa, marginals = margi,
#'  trawl.function = trawl.function)
rtrawl <- function(alpha, beta, times, marginal=rgamma, trawl.function=NA, trawl_fs=NA, trawl_fs_prim=NA, n, kappa = 0,
                    transformation=F, offset_shape=NULL, offset_scale=NULL, deep_cols=30){
  if(!is.na(trawl.function)){
    if(trawl.function %in% c("exp")){
      trawl_fs <- CollectionTrawl(times = times,
                                  params = list("rho"=rho), type = trawl.function)
      trawl_fs_prim <- CollectionTrawl(times = times,
                                       params = list("rho"=rho), type = trawl.function,
                                       prim = T)
    }else{
      stop(paste("trawl.function", trawl.function, "not yet implemented."))
    }

  }else{
    if(length(trawl_fs) != length(times)){
      stop('Wrong number of trawl functions compared to timestamps.')
    }
    if(length(trawl_fs_prim) != length(times)){
      stop('Wrong number of trawl primitives compared to timestamps.')
    }
  }

  if(transformation & (is.null(offset_scale) | is.null(offset_shape))){
    stop('When using marginal trf, indicate shape and scale offsets.')
  }

  if(!transformation){
    results <- TrawlSliceReconstruct(alpha = alpha,
                                      beta = beta+kappa,
                                      times = times,
                                      trawl_fs = trawl_fs,
                                      trawl_fs_prim = trawl_fs_prim,
                                      n = n, deep_cols)
  }else{
    results <- TrawlSliceReconstruct(alpha = offset_shape,
                                      beta = offset_scale,
                                      times = times,
                                      trawl_fs = trawl_fs,
                                      trawl_fs_prim = trawl_fs_prim,
                                      n = n, deep_cols)
  }

  return(results)
}

#' Simulation of extreme value path using trawl process. Transformed marginals
#' have scale parameter \code{1+kappa}.
#'
#' @param alpha Shape parameter.
#' @param beta Scale parameter.
#' @param kappa Additive constant to scale parameter \code{beta}.
#' @param times Vectors of discret times.
#' @param trawl_fs collection of trawl functions indexed on \code{times}.
#' @param trawl_fs_prim collection of trawl functions primitives indexed on
#'   \code{times}.
#' @param n Number of simulations (so far, only \code{n=1} is implemented).
#' @param transformation Boolean to apply marginal transform method. Default is
#'   False (F).
#' @param n_moments Number of finite moments for transformed marginals.
#' @param deep_cols Depth of reconstruction (columns). Default is 30.
#'
#' @returns Simulated path (size the same as times) of latent-trawl extreme
#'   value process.
#' @example TODO
rlexceed <- function(alpha, beta, kappa, times, trawl_fs, trawl_fs_prim, n, transformation, n_moments = 4, deep_cols=30){
  # TODO n is not implemented yet
  offset_shape <- n_moments+1
  offset_scale <- trf_find_offset_scale(alpha = alpha,
                                        beta = beta,
                                        kappa = kappa,
                                        offset_shape = n_moments+1)
  #print(offset_scale)
  # Generate Gamma
  gen_trawl <- rltrawl(alpha = alpha,
                       beta = beta,
                       kappa = 0.0,
                       times = times,
                       trawl_fs = trawl_fs,
                       trawl_fs_prim = trawl_fs_prim,
                       n = n,
                       transformation = transformation,
                       offset_shape = offset_shape,
                       offset_scale = 1)

  #return(gen_trawl)
  # Uniform threshold
  unif_samples <- runif(n=length(times)-deep_cols)
  if(n == 1){
    gen_exceedances <- rep(0, length(times)-deep_cols)
  }else{
    gen_exceedances <- matrix(0, nrow = length(times)-deep_cols, ncol = n)
  }

  #print(gen_trawl)
  prob_zero <- 1.0-exp(-kappa * gen_trawl)
  which_zero <- which(prob_zero >= unif_samples)
  if(transformation){
    gen_exceedances[-which_zero] <-  vapply(rexp(n = length(gen_trawl[-which_zero]), rate = gen_trawl[-which_zero]),
                                            FUN.VALUE = 1.0,
                                            FUN = function(x){return(trf_g(x, alpha = alpha,
                                                                           beta = beta,
                                                                           kappa = kappa,
                                                                           offset_scale = offset_scale,
                                                                           offset_shape = offset_shape))})
  }else{
    gen_exceedances[-which_zero] <-  rexp(n = length(gen_trawl[-which_zero]), rate = gen_trawl[-which_zero])
  }
  #mean(gen_exceedances)
  return(gen_exceedances)
}

# Various tests
# {
#   # Example
#   n_sims <- 50
#   times <- 1:400
#   kappa <- 0.3
#   alpha <- 3
#   beta <- 6
#   rho <- 0.3
#   n_moments <- 4
#
#   ## Find offset scale
#   offset_shape <- n_moments + 1
#   kappa / ((1+kappa/beta)^{alpha/offset_shape} - 1)
#   trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
#   offset_scale  <- trf_find_offset_scale(alpha = alpha, beta = beta, kappa = kappa, offset_shape = offset_shape)
#
#   cat("Prob non zero for non-trf",(1+kappa/beta)^{-alpha}, "\n")
#   cat("Prob non zero for trf",(1+kappa/offset_scale)^(-offset_shape), "\n")
#
#   ## Trawl process simulation
#   library(gPdtest)
#
#   ### Generating the functions
#   trawl_1 <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp", prim = F)
#   trawl_1_prim <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp", prim = T)
#
#   trl_slice <- trawl_slice_sets_not_optim(alpha = alpha, beta = beta, times = times, trawl_fs = trawl_1, trawl_fs_prim = trawl_1_prim, n = 1)
#
#   ### no transformation
#   gen_trawl <- rltrawl(alpha = alpha,
#                         beta = beta,
#                         times = times,
#                         n = 1,
#                         trawl_fs = trawl_1,
#                         trawl_fs_prim = trawl_1_prim,
#                         kappa = 0.0,
#                         transformation = F)
#   acf(gen_trawl, type = "covariance")
#   (alpha)/(beta)^2
#   hist(gen_trawl, probability = T)
#   lines(seq(0.01, 8, length.out = 200),dgamma(seq(0.01, 8, length.out = 200), shape = alpha, scale = 1/(beta)), col = "red")
#
#   par(mfrow=c(1,2))
#   #### ACF
#   acf(gen_trawl, main = paste("ACF trawl with rho =", rho))
#   lines(0:20, exp(-rho*0:20), col = "red")
#
#   #### distribution
#   plot(density(gen_trawl), "Marginal density of trawls")
#   lines(density(rgamma(n = 1000, shape = alpha, rate = beta)), col="red")
#
#   ### no transformation
#   (1+kappa/beta)^{-alpha}
#   1/(1+kappa)
#   par_ests_sims_no_trf <- matrix(0, ncol = 2, nrow = n_sims)
#   for(i in 1:n_sims){
#     gen_exc <- rlexceed(alpha = alpha,
#                         beta = beta,
#                         kappa = kappa,
#                         trawl_fs = trawl_1,
#                         trawl_fs_prim = trawl_1_prim,
#                         times = times,
#                         n = 1,
#                         transformation = F)
#     print(mean(gen_exc))
#     par_ests_sims_no_trf[i,] <- fExtremes::gpdFit(gen_exc, u =1e-6)@fit$par.ests
#   }
#   cat("mean:", (1+kappa/beta)^{-alpha}*(beta+kappa)/(alpha-1), "\n")
#
#   #### xi
#   1/alpha
#   boxplot(par_ests_sims_no_trf[,1])
#   abline(h=1/alpha, col = "red")
#   mean(par_ests_sims_no_trf[,1])
#   sd(par_ests_sims_no_trf[,1])
#
#   #### sigma
#   (beta+kappa)/alpha
#   boxplot(par_ests_sims_no_trf[,2])
#   abline(h=(beta+kappa)/alpha, col = "red")
#   mean(par_ests_sims_no_trf[,2])
#   sd(par_ests_sims_no_trf[,2])
#   # OBS: depending on fitting procedure used: over or under estimation happening gPdtest::gpd.fit and fExtremes::gpdFit
#
#
#   ### transformation
#   par_ests_sims_trf <- matrix(0, ncol = 2, nrow = n_sims)
#   for(i in 1:n_sims){
#     gen_exc_trf <- rlexceed(alpha = alpha,
#                             beta = beta,
#                             kappa = kappa,
#                             trawl_fs = trawl_1,
#                             trawl_fs_prim = trawl_1_prim,
#                             times = times,
#                             n = 1,
#                             transformation = T)
#     par_ests_sims_trf[i,] <- fExtremes::gpdFit(gen_exc_trf, u =1e-6)@fit$par.ests
#   }
#
#   #### xi
#   10
#   1/alpha
#   boxplot(par_ests_sims_trf[,1])
#   abline(h=1/alpha, col = "red")
#   mean(par_ests_sims_trf[,1])
#   sd(par_ests_sims_trf[,1])
#
#   #### sigma
#   (beta+kappa)/alpha
#   (beta)/alpha
#   boxplot(par_ests_sims_trf[,2])
#   abline(h=(beta)/alpha, col = "red")
#   mean(par_ests_sims_trf[,2])
#   sd(par_ests_sims_trf[,2])
# }
