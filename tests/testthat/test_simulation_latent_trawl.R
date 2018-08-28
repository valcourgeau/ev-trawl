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
  tt.prim <- TrawlExpPrimitive(t=times, rho=rho)

})


