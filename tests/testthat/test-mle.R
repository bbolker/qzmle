library(testthat)

set.seed(101)
d <- data.frame(x = 0:10, y = c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
if (requireNamespace("emdbook")) {
  data("ReedfrogPred", package="emdbook")
  rfp <- transform(ReedfrogPred, nsize = as.numeric(size), random = rnorm(48))
  form <- surv ~ dbinom(size = density, prob = exp(log_a) / (1 + exp(log_a) * h * density))
  form_logit <- surv ~ dbinom(size = density, prob = plogis(logit_a) / (1 + plogis(logit_a) * h * density))
  form_logit2 <- surv ~ dbinom(size = density,
                               prob = 1/(1 + exp(-logit_a)) / (1 + 1/(1+exp(-logit_a)) * h * density))
  form_logit3 <- surv ~ dbinom(size = density,
                               prob = 1/(1 + exp(-logit_a)) / (1 + 1/(1+exp(-logit_a)) * exp(log_h) * density))

  ## FIXME: does plogis() work too (or do we need to add it to deriv table?)

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
    suppressWarnings(qzfit <- qzmle::mle(form,
                 start = list(h = 4, log_a = 2),
                 parameters = list(log_a ~ poly(random)), data = rfp
                 ))
    ## qzfit2 <- qzmle::mle(form_logit,
    ## start = list(h = 4, logit_a = 2),
    ## parameters = list(logit_a ~ poly(random)), data = rfp
    ## )
    ## plogis() not in derivative table ...
    qzfit2 <- qzmle::mle(form_logit3,
                         start = list(log_h = log(4), logit_a = 2),
                         parameters = list(logit_a ~ poly(random)), data = rfp
                         )
    ## wants h -> 0 ... why???

    ## FIXME:: warnings about a>1; use logit link ???
    if (requireNamespace("bbmle")) {
      suppressWarnings(bbfit <- bbmle::mle2(form,
                           parameters = list(log_a ~ poly(random)),
                           start = list(log_a = c(2, 0), h = 4), data = rfp
                           ))
      expect_equal(round(bbmle::coef(bbfit), 1), round(qzfit$coefficients), 1)
      expect_equal(round(bbfit@min, 1), round(qzfit$minuslogl, 1))
    }
  })
} ## emdbook available



