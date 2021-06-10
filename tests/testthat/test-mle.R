library(testthat)

set.seed(101)
d <- data.frame(x = 0:10, y = c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
if (requireNamespace("emdbook")) {
  data("ReedfrogPred")
  rfp <- transform(emdbook::ReedfrogPred, nsize = as.numeric(size), random = rnorm(48))
  form <- surv ~ dbinom(size = density, prob = exp(log_a) / (1 + exp(log_a) * h * density))

  ## test if getting same result as bbmle

  test_that("Normal dist. works", {
    qzfit <- mle(y ~ dnorm(mean = ymean, sd = exp(logsd)),
                 start = list(ymean = mean(d$y), logsd = log(sd(d$y))),
                 data = d
                 )

    if (requireNamespace("bbmle")) {
      bbfit <- bbmle::mle2(y ~ dnorm(mean = ymean, sd = exp(logsd)),
                           start = list(ymean = mean(d$y), logsd = log(sd(d$y))),
                           data = d
                           )
      expect_equal(bbmle::coef(bbfit), qzfit$coefficients)
    }
  })

  test_that("Binomial dist. with Submodel works", {
    qzfit <- mle(form,
                 start = list(h = 4, log_a = 2),
                 parameters = list(log_a ~ poly(random)), data = rfp
                 )

    if (requireNamespace("bbmle")) {
      bbfit <- bbmle::mle2(form,
                           parameters = list(log_a ~ poly(random)),
                           start = list(log_a = c(2, 0), h = 4), data = rfp
                           )
      expect_equal(round(bbmle::coef(bbfit), 1), round(qzfit$coefficients), 1)
      expect_equal(round(bbfit@min, 1), round(qzfit$minuslogl, 1))
    }
  })
} ## emdbook available
