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
    # qzfit2 <- qzmle::mle(form_logit3,
    #                      start = list(log_h = log(4), logit_a = 2),
    #                      parameters = list(logit_a ~ poly(random)), data = rfp
    #                      )
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


set.seed(1001)
lymax <- c(0,2)
lhalf <- 0
x <- sort(runif(200))
g <- factor(sample(c("a","b"),200,replace=TRUE))
y <- rnbinom(200,mu=exp(lymax[g])/(1+x/exp(lhalf)),size=2)
d2 <- data.frame(x,g,y)

test_that("Negative Binomial works", {

  fit3qz <- qzmle::mle(y~dnbinom(mu=exp(lymax)/(1+x/exp(lhalf)),
                                 size=exp(logk)),
                       parameters=list(lymax~g), data=d2,
                       start=list(lymax=0,lhalf=0,logk=0))

  if (requireNamespace("bbmle")) {
    fit3bb <- bbmle::mle2(y~dnbinom(mu=exp(lymax)/(1+x/exp(lhalf)),size=exp(logk)),
                   parameters=list(lymax~g),data=d2,
                   start=list(lymax=0,lhalf=0,logk=0))
  }
  expect_equal(bbmle::coef(fit3bb), coef(fit3qz), tolerance = 1e-6)
  }
  )

d <- data.frame(x = 0:10,
               y = c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
attach(d)

LL <- function(x, y, ymax=15, xhalf=6, log=FALSE) {
  return(-sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE)))
}
test_that("Cannot use data in local environment", {
  expect_error(qzmle::mle(LL), "missing `data` argument")
})

test_that("Poisson works", {
  fit0qz <- qzmle::mle(y~dpois(lambda=ymean),
                     start=list(ymean=mean(d$y)),
                     data=d)
  if (requireNamespace("bbmle")) {
    fit0bb <- bbmle::mle2(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
  }
  expect_equal(coef(fit0qz), bbmle::coef(fit0bb), tolerance = 1e-5)
  expect_equal(vcov(fit0qz), bbmle::vcov(fit0bb), tolerance = 1e-5)
  expect_equal(logLik(fit0qz), bbmle::logLik(fit0bb), tolerance = 1e-4)
})


if (requireNamespace("bbmle")) {
  fit1bb <- bbmle::mle2(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
      start=list(lymax=0,lhalf=0),
      data=d,
      parameters=list(lymax~1,lhalf~1))

  test_that("Poisson with more than on parameter works", {
  fit1qz <- qzmle::mle(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
      start=list(lymax=0,lhalf=0),
      data=d,
      parameters=list(lymax~1,lhalf~1))

  expect_equal(unname(bbmle::coef(fit1bb)), unname(coef(fit1qz)), tolerance = 1e-6)
  })

  test_that("(TMB) Poisson with more than on parameter works", {
  fit1qz_tmb <- qzmle::mle(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
      start=list(lymax=0,lhalf=0),
      data=d,
      parameters=list(lymax~1,lhalf~1),
      method = "TMB")
  expect_equal(unname(bbmle::coef(fit1bb)), unname(coef(fit1qz_tmb)), tolerance = 1e-6)
  unlink("template.*")
  })
}
