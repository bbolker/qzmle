## let's say I want to fit a functional response model

killed ~ dbinom(size=initial_N,
                prob=exp(log_a)*initial_N/(1+exp(log_a)*h*initial_N))
## but what if I want log_a to be a quadratic function of predator size?
parameters=list(log_a ~ poly(predsize,2))
parameters=list(log_a ~ 1 + predsize + I(predsize^2))
parameters=list(log_a ~ predsize)  ## automatically add an intercept
## internally we need to set up a model matrix ...
X <- model.matrix(~poly(predsize,2), data=data)
## log_a coefficients will now be a _vecotr_ of three values
log_a <- X %*% log_a_coeffs


## when you create the LL function, you need to set up a list of X matrices
## for all of the parameters included in the 'parameters=' argument
## when you evaluate the LL function, you need to pull out the sub-lists of parameters
## for each component and %*% them by the corresponding X to get the value that's
## used in evaluating the distribution
## ... AND ... you need another chain-rule step!
## D(NLL, { log.a.intercept, log.a.slope, log.a.quad, h }) =
## D(dbinom/prob)*D(prob/(log_a))* {D(log_a)/(log.a.intercept),
##                                  D(log_a)/(log.a.slope),
##                                  D(log_a)/(log.a.quad)}
## {D(dnorm/mean), D(dnorm/sd)}


parameter_parse <- function(parameters, data){
  ## assume only one submodel for now
  formula <- parameters[[1]]
  RHS <- formula[-2] ## "~ size"

  ## set up model
  X <- model.matrix(RHS, data=data)

  ## convert X into a list
  val_list <- split(X, rep(1:ncol(X), each = nrow(X)))
  names(val_list) <- colnames(X)
}


library(bbmle)
library(emdbook)
mle2(surv ~ dbinom(size=density,
                  prob=1/(1+exp(log_a)*h*density)),
    parameters=list(log_a~size),
    start=list(log_a=0,h=1),
    data=ReedfrogPred)


model.matrix(~size, data=ReedfrogPred)


     ## relist() ???
orig_list <- list(a=c(1,2,3),b=4)
vec <- unlist(orig_list)
relist(vec, orig_list)

## orig arguments
list(a=c(intercept, slope, quad),
     h=2)
## mlefun must take a VECTOR of params
## because that's what optim() gives it
relist(p,...)


## *when* we do this in TMB we will want
## to pass some kind of vector of indices
## saying which elements of the param vec
## correspond to individual parameters in
## the model



