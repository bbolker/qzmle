## 1. link-transform response (probably not), e.g. log(y) ~ ...
## 2. inverse-link-transform prediction (as in GLMs), e.g.
##    y ~ exp(beta0+beta1*x)
## 3. inverse-link-transform *coefficients*. Example:
##  type-2 functional response
##   consumption ~  dbinom(N=exposed, prob = a/(1+a*h*N))
##   in order to ensure a, h > 0,  a = exp(log_a), h = exp(log_h)
##   links = list(a="log", h="log")
##
make.link("crazylink")
##  TMB code ...
##  PARAMETER(log_a);
##  PARAMETER(log_h);
##  Type a = exp(log_a);
##  Type h = exp(log_h);
##  nll = -sum(dbinom(exposed, a/(1+a*h*N), true));
##
## R code version:
##     -sum(dbinom(exposed,exp(log_a)/(1+exp(log_a)*exp(log_h)*N),
##           log=TRUE))
##  ALTERNATIVE: implement another chain rule step
##     d(prob)/d(log_a) = d(prob)/da * da/d(log_a)

## so if dvec = a gradient vector with response to the
##   response-scale prameters (dL/da, dL/dh)
## then each element of the vector needs to be multiplied by the
##  chain-rule element ...
grad <- c(1,2,1,1)
linkvec <- c("log","identity","logit","cloglog")
## names of all allowed links (except 'identity')
all_links <- c("log","logit","cloglog","sqrt","inverse","log10")
trans_parnames <- function(p) {
    regex <- sprintf("(%s)_", paste(all_links,collapse="|"))
    gsub(regex,"",p)
}
## params on link scale
orig_parvec <- c(log_a=-1, h = 2, logit_q=3, cloglog_r=-2)
parvec <- numeric(length(orig_parvec))
names(parvec) <- trans_parnames(names(orig_parvec))
g <- numeric(length(orig_parvec))
for (i in seq_along(linkvec)) {
    mm <- make.link(linkvec[i])  ## get the whole 'link' object
    parvec[i] <- mm$linkinv(orig_parvec[i])
    g[i] <- 1/mm$mu.eta(parvec[i])
}

grad <- grad*g ## apply link-function chain rule ...

deriv(expression(log(x)),"x")
D(expression(log(x)),"x")
D(expression(x), "log(x)")
library(Deriv)
Deriv(expression(log(x)),"x")
Deriv(expression(log(x)),exp("x"))
## D(f^{-1}(x)) = 1/D(f(x))
## D(exp(x)/x) = exp(x)
## D(log(x)/x) = 1/x = 1/exp(log(x)) = 1/D(exp(log(x)))/D(log(x))
##  = 1/( d(f^{-1}(y))/y ) = 1/mu.eta(y) { y is on the link scale }

## inverse-link for log = exp
## log_a = 2
## a = exp(2)
## d(a)/d(log_a)

## evaluate at the value of log_a
## 1/deriv(exp) = deriv(log)
## should be the same as the deriv d(log(y))/dy 
## of the LINK function at the RESPONSE value
## i.e. 1/y = 1/invlink(x) = 1/exp(x)
## x == log(a) [LINK scale] , y == a [RESPONSE scale]
## LOG: LINK y -> x,  EXP INVLINK x -> y

log_a <- 2
## d(eta)/d(mu) = d(log_a)/da
1/make.link("log")$mu.eta(log_a)
## d(eta)/d(mu) = d(log(y))/dy = 1/y
1/(make.link("log")$linkinv(log_a))

logit_a <- 2
1/make.link("logit")$mu.eta(logit_a)
## derivative

##
##  derivation (ignore for now)
##   where a = attack rate, h = handling time
##   time taken to find a new prey item = a*N
##   time taken to handle prey once you find it = h
##   (handling time = pursuit + killing + digestion + ...)
##   total time between attacks = a*N + h
##   RATE of acquiring new prey = 1/(a*N+h)

##
b <- body(make.link)
b[[2]] ## guts of the function
linknames <- names(b[[2]])
linknames <- linknames[nzchar(linknames)]
