library(testthat)

set.seed(101)
d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))

## test if getting same result as bbmle

test_that("Normal dist. works", {
  qzfit <- mle(y~dnorm(mean=ymean, sd=exp(logsd)),
             start=list(ymean=mean(d$y), logsd=log(sd(d$y))),
             data=d)

  if (requireNamespace("bbmle")) {
  bbfit <- bbmle::mle2(y~dnorm(mean=ymean, sd=exp(logsd)),
                               start=list(ymean=mean(d$y), logsd=log(sd(d$y))),
                               data=d)
  expect_equal(bbmle::coef(bbfit), qzfit$coefficients)
  }
})


