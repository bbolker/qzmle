library(testthat)

## set up general purpose objects for testing
set.seed(101)
dd <- data.frame(
  y = rpois(100, lambda = 2),
  latitude = rnorm(100)
)
## manual function
## pois needs x in its environment
x <- dd$y
pois <- function(lambda) {
  x * log(lambda) - lambda - lfactorial(x)
}
normal <- function(mean, sd) {
  -log((2 * pi)^0.5) - log(sd) - (x - mean)^2 / (2 * sd^2)
}
## MUST use this (rather than dpois) for now:
## FIXME: can we get numDeriv() to work with a multi-parameter model?

test_that("Poisson dist. with single parameter lambda works", {
  ## formula and data
  form <- y ~ dpois(lambda = b0 * latitude^2)

  ## my function
  ff <- mkfun(form, start = list(b0 = 3), data = dd)
  expect_equal(names(ff), c("start", "fn", "gr"))

  ## manually setting parameters
  x <- dd$y ## pois looks for x in its environment
  lambda <- function(b0) {
    b0 * dd$latitude^2
  }

  expect_equal(
    ff$fn(c(b0 = 3)),
    -sum(pois(lambda(b0 = 3)))
  )

  b0 <- 3
  lambda <- b0 * dd$latitude^2
  if (requireNamespace("numDeriv")) {
    expect_equal(
      ff$gr(c(b0 = 3))[["lambda.b0"]],
      -sum(numDeriv::grad(pois, lambda) * dd$latitude^2)
    )
  }
})

test_that("Poisson with two parameters", {
  ## formula and data
  form2 <- y ~ dpois(lambda = b0 + b1 * latitude^2)

  ## my function
  ff2 <- mkfun(form2, start = list(b0 = 3, b1 = 2), data = dd)

  ## manually setting parameters
  lambda2 <- function(b0, b1) {
    b0 + b1 * dd$latitude^2
  }
  L <- lambda2(b0 = 3, b1 = 2)
  expect_equal(
    -sum(pois(L)),
    -sum(dpois(dd$y, L, log = TRUE))
  )

  expect_equal(
    ff2$fn(c(b0 = 3, b1 = 2)),
    -sum(pois(L))
  )

  if (requireNamespace("numDeriv")) {
    ## compute gradient by hand: have to be a little clever to
    ##  get scalar results in both L and x
    fullgrad <- rep(NA, length(L))
    for (i in seq_along(fullgrad)) {
      assign("x", dd$y[i], environment(pois))
      fullgrad[i] <- numDeriv::grad(pois, L[i])
    }
    assign("x", dd$y, environment(pois))

    ## chain rule by hand: dL/db0=1, dL/db1=dd$latitude^2
    ## d(loglik)/db0 = d(loglik)/d(lambda)*d(lambda)/d(b0)
    expect_equal(
      unname(ff2$gr(c(b0 = 3, b1 = 2))),
      c(-sum(1 * fullgrad), -sum(dd$latitude^2 * fullgrad))
    )
  }
})

test_that("normal with single parameter mean", {
  ## formula and data
  form3 <- y ~ dnorm(mean = b0 * latitude^2, sd = 1)
  ff3 <- mkfun(form3, start = list(b0 = 3), data = dd)

  m <- function(b0) {
    b0 * dd$latitude^2
  }
  expect_equal(
    ff3$fn(c(b0 = 3)),
    -sum(normal(m(b0 = 3), sd = 1))
  )
})

test_that("normal with two parameter mean", {
  ## formula and data
  form4 <- y ~ dnorm(mean = b0 + b1 * latitude^2, sd = 1)
  ff4 <- mkfun(form4, start = list(b0 = 1, b1 = 2), data = dd)
  ff4$fn(c(b0 = 1, b1 = 2))

  m <- function(b0, b1) {
    b0 + b1 * dd$latitude^2
  }
  M <- m(b0 = 1, b1 = 2)
  expect_equal(
    ff4$fn(c(b0 = 1, b1 = 2)),
    -sum(normal(m(b0 = 1, b1 = 2), sd = 1))
  )
})


## optim(par=<starting parameters>, fn= ...$fn, gr = ... $gr, method="BFGS")
## default method is Nelder-Mead, which ignores gradient
## store the reference value so you can test against it
## dput() prints out a full

fr <- function(x) { ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
  x1 <- x[1]
  x2 <- x[2]
  c(
    -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
    200 * (x2 - x1 * x1)
  )
}
xx <- optim(c(-1.2, 1), fr)
## dput(xx)
ref_val <- list(
  par = c(1.00026013872567, 1.00050599930377), value = 8.82524109672275e-08,
  counts = c(`function` = 195L, gradient = NA), convergence = 0L,
  message = NULL
)
expect_equal(xx, ref_val)
