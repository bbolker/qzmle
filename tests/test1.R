library(devtools)
use_r("tmb")
load_all()
check()


## Test 1 - 1 model param, distribution with 1 param -------
form <- y ~ dpois(lambda = b0 * latitude^2)
set.seed(101)
dd <- data.frame(
  y = rpois(100, lambda = 2),
  latitude = rnorm(100)
)

## debug(mkfun)
ff <- mkfun(form, data=dd)
## debug(ff)
ff$fn(c(b0=1)) ##417.6891
ff$gr(c(b0=1)) ## -113.5597

##check fn() with finite differences
lambda <- function(b0) {b0*dd$latitude^2}
x = dd$y
-sum(x*log(lambda(b0=1))-lambda(b0=1)-lfactorial(x)) ## 417.6891

##check gr() with numDeriv grad()
b2=1
x2=dd$y
lambda2 = b2*dd$latitude^2
pois <- function(lambda2){x2*log(lambda2)-lambda2-lfactorial(x2)}
-sum(grad(pois,lambda2)*dd$latitude^2) ##-113.5597






## PDF param <- multiple model parameters
form2 <- y ~ dpois(lambda = b0 + b1 * latitude^2)
## ff <- mkfun(form2)
## what happens if we forget data? error ...
ff2 <- mkfun(form2, data=dd)
ff2$fn(c(b0=1,b1=2))
undebug(ff$gr)
ff2$gr(c(b0=1,b1=2))

## testing
ff$gr(c(b1=2)) == ff2$gr(c(b0=0, b1=2))[["b1"]]

-sum(dpois(dd$y, 1 + 2*dd$latitude^2, log=TRUE))



### more implementations ---------------
debug(ff$gr)
ff$gr(c(b0=1,b1=2))


form3 <- y ~ dnorm(mean = b0 * latitude^2, sd = sd)
form4 <- y ~ dnorm(mean = b0 + b1 * latitude^2, sd = sd)




##### JUNK -------------

-sum(dpois(dd$y, lambda = 1 * dd$latitude^2, log = TRUE))

## deriv of -neg log likelihood of
-sum(dpois(y, lambda = b0 * latitude^2, log = TRUE))

## -sum(deriv(loglik(lambda,x))/lambda *
##     deriv(lambda/parameters)


##
deriv(lambda_nll, "lambda")
library(Deriv)

eval(d1, list(x = 1, lambda = 2))

