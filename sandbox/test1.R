## library(devtools)
## load_all()
## check()


## form -> objective function (i.e. a function that computes -sum(dpois(y, lambda=b0*latitude^2))
##  where the variables in the formula are *evaluated* in an environment that
##  includes the current parameter
##  values and information stored in a 'data' variable


## Test 1 - 1 model param, dist with 1 param -------
form <- y ~ dpois(lambda = b0 * latitude^2)
set.seed(101)
dd <- data.frame(
  y = rpois(100, lambda = 2),
  latitude = rnorm(100)
)

set.seed(101)
dd2 <- data.frame(
  y = rpois(100, lambda = 2),
  latitude = rnorm(100),
  sd=1
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
##----------
-sum(dpois(dd$y, 1*dd$latitude^2, log=TRUE)) == ff$fn(c(b0=1))
## TRUE

##check gr() with numDeriv grad()
b0=1
lambda = b0*dd$latitude^2
pois <- function(lambda){x*log(lambda)-lambda-lfactorial(x)}
### d/d(lambda) * d lambda / d(b2)
-sum(numDeriv::grad(pois,lambda)*dd$latitude^2) ## -113.5597



## Test 2 - multiple model parameters, dist with 1 param -------------
form2 <- y ~ dpois(lambda = b0 + b1 * latitude^2)
ff <- mkfun(form2) ## what happens if we forget data? output error message ...
ff2 <- mkfun(form2, data=dd)
ff2$fn(c(b0=1,b1=2)) ## 215.3971
ff2$gr(c(b0=1,b1=2)) ## b0 = -4.823218 b1=43.851884

## testing
ff$gr(c(b0=2)) == ff2$gr(c(b0=0, b1=2))[["b1"]] ## TRUE; -8.559725

###Check fn() with finite diff
-sum(dpois(dd$y, 1 + 2*dd$latitude^2, log=TRUE))  == ff2$fn(c(b0=1,b1=2))
## TRUE

##Check gr() with numDeriv grad()
b0=1
b1=2
lambda = b0 + b1*dd$latitude^2
pois <- function(lambda){x*log(lambda)-lambda-lfactorial(x)}
### d/d(lambda) * d lambda / d(b0)
-sum(numDeriv::grad(pois,lambda)*1) #-4.823218
ff2$gr(c(b0=1,b1=2))[['b0']] #-4.823218

### d/d(lambda) * d lambda / d(b1)
-sum(numDeriv::grad(pois,lambda)*dd$latitude^2) #43.85188
ff2$gr(c(b0=1,b1=2))[['b1']] #43.85188


### Test 3 --------------
form4 <- y ~ dnorm(mean = b0 + b1 * latitude^2, sd = 1)
ff4 <- mkfun(form4, data=dd)
ff4$fn(c(b0=1, b1=2)) == -sum(dnorm(dd$y, mean=1 + 2*dd$latitude^2, sd=1, log=TRUE))
## TRUE





### more implementations ---------------
form3 <- y ~ dnorm(mean = b0 * latitude^2, sd = sd)

form4 <- y ~ dnorm(mean = b0 + b1 * latitude^2, sd = sd)



-sum(dnorm(dd$y, mean=1 + 2*dd$latitude^2, sd=1, log=TRUE))

##### JUNK -------------

-sum(dpois(dd$y, lambda = 1 * dd$latitude^2, log = TRUE))

## deriv of -neg log likelihood of
-sum(dpois(y, lambda = b0 * latitude^2, log = TRUE))

## -sum(deriv(loglik(lambda,x))/lambda *
##     deriv(lambda/parameters)

## confirming dnorm loglik func
mean1=1 + 2*dd$latitude^2
sd=1
-sum(dnorm(dd$y, mean=1 + 2*dd$latitude^2, sd=1, log=TRUE)) == -sum(-log(2*pi)/2-log(sd) - (dd$y-mean1)^2/(2*sd^2))



##
deriv(lambda_nll, "lambda")
library(Deriv)

eval(d1, list(x = 1, lambda = 2))

