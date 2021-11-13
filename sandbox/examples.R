## examples to make work with qzmle
library(qzmle)
library(bbmle)
library(testthat)

## examples from ?bbmle::mle2

d <- data.frame(x = 0:10,
               y = c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
attach(d) ## !!!!
## in general it is best practice to use the `data' argument,
##  but variables can also be drawn from the global environment
## LL <- function(ymax=15, xhalf=6)
## -sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE))
LL <- function(x, y, ymax=15, xhalf=6, log=FALSE) {
  return(-sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE)))
}
add_logl("LL", expression(-sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE))),
         c("ymax", "xhalf"))
## uses default parameters of LL
## data from global environment
fit <- mle2(LL)
expect_error(qzmle::mle(LL), "function does not use data")

## TASK 1. (a) Document the fact that qzmle::mle() doesn't work
## with variables taken from the global environment
## (b) Add an issue to possibly implement this ...
## There's a legitimate argument for **not** implementing this.
## (Because parameters and data are grouped before passing to optim)

## TASK 2.
## (a) Document that qzmle::mle() doesn't work with user-specified NLL functions.
## (b) Add an issue to implement this.
fit <- mle2(LL, data=d)
fitqz <- qzmle::mle(LL, data = d)

fit0 <- mle2(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
fit0qz <- qzmle::mle(y~dpois(lambda=ymean),
                     start=list(ymean=mean(d$y)),
                     data=d)
## Cant put minuslogl inside bc minuslogl = opt$value
##optimHess(coef(fit0qz), fn = fit0qz$minuslogl)
expect_equal(coef(fit0), coef(fit0qz))
expect_equal(vcov(fit0), vcov(fit0qz), tol=1e-5)
## TASK: return the negative log-likelihood _function_ as one component of a fitted mle object.


## TASK: improve log likelihood method
## FROM BBMLE:

## function (object, ...) {
##     if (length(list(...)))
##         warning("extra arguments discarded")
##     ## get negative log-likelihood
##     val <- -object@min
##     ## get df/number of coefficients
##     attr(val, "df") <- length(object@coef)
##     ## get number of observations
##     attr(val, "nobs") <- attr(object, "nobs")
##     class(val) <- "logLik"
##     val
## }

## don't print, just return an appropriately structured object of
## class "logLik" ...

## anova(fit0,fit)
summary(fit0)
logLik(fit0)
vcov(fit0)

summary(fit0qz)

expect_equal(logLik(fit0qz), logLik(fit0), tolerance = 1e-4)

vcov(fit0qz)
## Able to assign to variable
temp <- vcov(fit0qz)
temp <- logLik(fit0qz)


fit1qz <- qzmle::mle(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
      start=list(lymax=0,lhalf=0),
      data=d,
      parameters=list(lymax~1,lhalf~1))

fit1qz_tmb <- qzmle::mle(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
      start=list(lymax=0,lhalf=0),
      data=d,
      parameters=list(lymax~1, lhalf~1),
      method="TMB")

fit1bb <- mle2(y~dpois(lambda=exp(lymax)/(1+x/exp(lhalf))),
      start=list(lymax=0,lhalf=0),
      data=d,
      parameters=list(lymax~1,lhalf~1))

## TMB parameters are named "_param"
all.equal(unname(coef(fit1bb)),
          unname(coef(fit1qz_tmb)),
          tolerance = 1e-6)

set.seed(1001)
lymax <- c(0,2)
lhalf <- 0
x <- sort(runif(200))
g <- factor(sample(c("a","b"),200,replace=TRUE))
y <- rnbinom(200,mu=exp(lymax[g])/(1+x/exp(lhalf)),size=2)
d2 <- data.frame(x,g,y)

fit3bb <- mle2(y~dnbinom(mu=exp(lymax)/(1+x/exp(lhalf)),size=exp(logk)),
             parameters=list(lymax~g),data=d2,
             start=list(lymax=0,lhalf=0,logk=0))

fit3qz <- qzmle::mle(y~dnbinom(mu=exp(lymax)/(1+x/exp(lhalf)),
                               size=exp(logk)),
                     parameters=list(lymax~g), data=d2,
                     start=list(lymax=0,lhalf=0,logk=0))

compare(coef(fit3bb), coef(fit3qz), tolerance = 1e-6)
all.equal(coef(fit3bb), coef(fit3qz), tol=1e-6)

## TASK 1: copy these examples over to test files,
## with expect_equal() instead of all.equal()
## (expect_error() where used above)

## TASK 2: try TMB with random effects!

## TASK 3: when parameters= is specified, parameter names should be
## coefficient names should match between mle and TMB versions
## parameter *vector* is still called lymax_params (inside TMB),
## but the names associated with the coefficients should be
## lymax.(Intercept) etc..

## for each model matrix, paste the parameter name together with the column names
##  of the model matrix
##  lymax . (Intercept)
## paste(param , colnames(X1), sep = ".")
