f <- function(a,b,c=2,d=3) {
    cat(a,b,c,d,"\n")
}
    
f(1,2,3,4)
f(0,0)
f(c=3,2,1)

y ~ dbinom(size=N,prob=a/(1+a*h*N))
dbinom(y,N,a/(1+a*h*N))

dbinom


library(bbmle)
devtools::load_all("..")
dd <- data.frame(y=rpois(100,lambda=1))

fun1 <- mkfun(y~dpois(exp(lambda)), data=dd)
fun2 <- mkfun(y~dnorm(mean = b0 + b1 * latitude^2, sd = 1), data=dd)
form <- surv ~ dbinom(size = density, prob = exp(log_a)/(1 + exp(log_a)*h*density))
rfp <- transform(emdbook::ReedfrogPred, nsize=as.numeric(size))
fun3 <- mkfun(form,parameters=list(log_a~poly(nsize)),data=rfp)
fun4 <- mkfun(form,parameters=list(log_a~poly(nsize), h~1),data=rfp)
debug(parameter_parse)

fun3$fn(c(log_a=0,h=1))
fun3$fn(c(log_a=1,h=1))

data("ReedfrogPred",package="emdbook")
m1 <- bbmle::mle2(form,parameters=list(log_a~poly(nsize,1)),
            data=rfp,
            start=list(log_a=0,h=1))
## the negative log-likelihood function is saved as @minuslogl in the fitted object
## the other trick is that we are automatically expanding the start value for a sub-model
## by taking the named argument as the first element, all other elements of the sub-model
## parameter vector are 0
m1@minuslogl(0,0,1)  ## log_a=0 -> (0,0)
m1@minuslogl(1,0,1)  ## log_a=1 -> (1,0)
formals(m1@minuslogl)

m1@minuslogl(1,2,1)

orig_pars <- list(log_a=c(NA,NA), h=NA)
pvec <- c(1,0,2)
relist(pvec,orig_pars)
##
## before defining $fn and $gr:
##    set up Xlist (i.e. calling model.matrix() for each piece)
##    make sure you have the orig_pars list so that you can use it to relist()
##     in step 1 *within* $fn and $gr
##     these object will live in the _environment_ of $fn and $gr, so they'll
##        be available
## internally within $fn (or $gr):
##  1. relist() the parameter vector to make it back into a structured list
##   i.e. relist(c(1,0,2), orig_pars) → parlist = list(log_a=c(1,0), h=2)
##  2. for $fn: loop over submodels (par = name of top-level parameter)
##      a. compute Xlist[[par]] %*% parlist[[par]]
##          Xlist = model matrix for submodel for parameter 'par'
##          parlist is the 'restructured' list from (a)
##      b. substitute these values back into the 'top-level' parameter list
##      c. evaluate the original expression with these values
##
##     for $gr:
##       
##  for each submodel
##  find the 

## environment example
## 1. working with global environment
f <- function() {
    print(a+1)
}
try(f())
a <- 5
f()
environment(f)

g <- function() {
    a <- 107
    tmpf <- function() {
        print(a+1)
    }
    return(tmpf)
}

f2 <- g()
f2()
ls(environment(f2))
ls(parent.env(environment(f2))) ##
parent.env(environment(f2))
