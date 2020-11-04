#library(numDeriv)

test_that("Poisson dist. with single parameter lambda works", {
  #formula and data
  form <- y ~ dpois(lambda = b0 * latitude^2)
  dd <- data.frame(y = rpois(100, lambda = 2),
    latitude = rnorm(100))

  # my function
  ff <- mkfun(form, data=dd)

  # manually setting parameters
  pois <- function(lambda){x*log(lambda)-lambda-lfactorial(x)}
  lambda <- function(b0) {b0*dd$latitude^2}
  x <- dd$y
  expect_equal(ff$fn(c(b0=3)),
               -sum(pois(lambda(b0=3))))

  b0=3
  lambda = b0*dd$latitude^2
  if (requireNamespace("numDeriv")) {
      expect_equal(ff$gr(c(b0=3))[['b0']],-sum(numDeriv::grad(pois,lambda)*dd$latitude^2))
  }
})

