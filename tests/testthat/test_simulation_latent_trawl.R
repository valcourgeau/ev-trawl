# source("./R/simulation_latent_trawl.R")

context("Trawl functions")

test_that("Exp trawl function rho positive", {
  t <- 10
  rho <- -0.2
  expect_error(TrawlExp(t = t, rho = rho))
})

test_that("Exponential trawl function values", {
  t <- 10
  rho <- 0.2

  t.f.exp <- TrawlExp(t = t, rho = rho)
  expect_equal(t.f.exp(t), 1)
  expect_equal(t.f.exp(-Inf), 0)
  expect_equal(t.f.exp(t + 10), 0)
})

test_that("Exp primitive trawl function rho positive", {
  t <- 10
  rho <- -0.2
  expect_error(TrawlExpPrimitive(t = t, rho = rho))
})

test_that("Exponential primitive trawl function values", {
  t <- 10
  rho <- 0.2

  t.f.exp <- TrawlExpPrimitive(t = t, rho = rho)
  expect_equal(t.f.exp(t), 1/rho)
  expect_equal(t.f.exp(-Inf), 0)
  expect_equal(t.f.exp(t + 10), 1/rho)
})

test_that("CollectionTrawl output length", {
  times <- 1:500/500
  rho <- 0.2
  expect_equal(length(CollectionTrawl(times, params = list("rho"=rho), type = "exp")),
               length(times))
})

test_that("CollectionTrawl output exp trawl function", {
  times <- 10
  rho <- 0.2
  collection.trawl <- CollectionTrawl(times, params = list("rho"=rho), type = "exp")
  exp.t.function <- collection.trawl[[1]]
  tt.fct <- TrawlExp(t=times, rho=rho)
  tt.xs <- seq(0, 20, by=0.5)

  tt.ys <- tt.fct(tt.xs)
  tt.zs <- exp.t.function(tt.xs)

  expect_equal(tt.ys, tt.zs)
})

context("Slice area and reconstruction")
test_that("Slice area value", {
  timestamps <- 10:12
  rho <- 0.2
  tt.prim <- TrawlExpPrimitive(t=timestamps[1], rho=rho)
  theoretical.val <- 1/rho - exp(rho*(timestamps[1]-timestamps[2]))/rho
  expect_equal(SliceArea(i=1, j=1, times = timestamps, trawl_f_prim = tt.prim), theoretical.val)
})

test_that("Reconstruction with wrong marginals", {
  set.seed(42)

  # params:
  alpha <- 3
  beta <- 2
  kappa <- 0
  n.timestamps <- 10
  times <- 1:n.timestamps
  n <- 1

  rho <- 0.2
  trawl.fs <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp")
  trawl.fs.prim <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp", prim = T)

  # marginals:
  marg <- "wrong marg"

  expect_error(TrawlSliceReconstruct(alpha = alpha, beta = beta, times = times,
                                  marg.dist = marg, n = n, trawl_fs = trawl.fs,
                                  trawl_fs_prim = trawl.fs.prim))
})

test_that("Reconstruction using SliceArea for Gamma marginals", {
  set.seed(42)

  # params:
  alpha <- 3
  beta <- 2
  kappa <- 0
  n.timestamps <- 3000
  times <- 1:n.timestamps
  n <- 1

  rho <- 0.2
  trawl.fs <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp")
  trawl.fs.prim <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp", prim = T)

  # marginals:
  marg <- "gamma"

  values <- TrawlSliceReconstruct(alpha = alpha, beta = beta, times = times,
                        marg.dist = marg, n = n, trawl_fs = trawl.fs,
                        trawl_fs_prim = trawl.fs.prim)
  to.compare <- rgamma(shape=alpha, rate=beta, n = n.timestamps)
  expect_gte(ks.test(values, "pgamma", alpha, beta)$p.value, 0.05)
})

test_that("Reconstruction using SliceArea for Gaussian marginals", {
  set.seed(42)

  # params:
  alpha <- 3
  beta <- 2
  kappa <- 0
  n.timestamps <- 3000
  times <- 1:n.timestamps
  n <- 1

  rho <- 0.2
  trawl.fs <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp")
  trawl.fs.prim <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp", prim = T)

  # marginals:
  marg <- "gaussian"

  values <- TrawlSliceReconstruct(alpha = alpha, beta = beta, times = times,
                                  marg.dist = marg, n = n, trawl_fs = trawl.fs,
                                  trawl_fs_prim = trawl.fs.prim)
  to.compare <- rnorm(n = n.timestamps, alpha, beta)
  expect_gte(ks.test(values,to.compare)$p.value, 0.05)

  set.seed(42)
  marg <- "normal"

  values <- TrawlSliceReconstruct(alpha = alpha, beta = beta, times = times,
                                  marg.dist = marg, n = n, trawl_fs = trawl.fs,
                                  trawl_fs_prim = trawl.fs.prim)
  to.compare <- rnorm(n = n.timestamps, alpha, beta)
  expect_gte(ks.test(values,to.compare)$p.value, 0.05)
})

test_that("Reconstruction using SliceArea for GIG marginals", {
  set.seed(42)

  # params:
  alpha <- 3
  beta <- 2
  kappa <- 0
  n.timestamps <- 3000
  times <- 1:n.timestamps
  n <- 1

  rho <- 0.2
  trawl.fs <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp")
  trawl.fs.prim <- CollectionTrawl(times = times, params = list(rho=rho), type = "exp", prim = T)

  # marginals:
  marg <- "gig"

  values <- TrawlSliceReconstruct(alpha = alpha, beta = beta, times = times,
                                  marg.dist = marg, n = n, trawl_fs = trawl.fs,
                                  trawl_fs_prim = trawl.fs.prim)
  to.compare <- rgig(n = n.timestamps, lambda = 0, chi = alpha, psi = beta)
  expect_gte(ks.test(values, to.compare)$p.value, 0.05)
})



