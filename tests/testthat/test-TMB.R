library(testthat)

set.seed(123)
x <- runif(20, 1, 10)
y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
form <- y ~ dnorm(b0 + b1 * x, log_sigma)
links <- c(b0 = "identity", b1 = "identity", sigma = "log")
start <- list(b0 = 1, b1 = 2, log_sigma = sd(y))

test_that("data passing to template", {
  data_list <- TMB_template(form, start = start,
                       links = links, data = list(x = x, y = y))

  expect_output(str(data_list), "List of 2")

})
