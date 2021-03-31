library(TMB)
compile("./submodel.cpp")
dyn.load(dynlib("./submodel"))

set.seed(101)
rfp <- transform(emdbook::ReedfrogPred,
                 nsize=as.numeric(size), random=rnorm(48))

form <- surv ~ dbinom(size = density,
                      prob = exp(log_a)/(1 + exp(log_a)*h*density))

##bbmle
bbmle::mle2(form,start=list(log_a=c(2,1), h=4),
            parameters=list(log_a~poly(nsize)),data=rfp)

## qzmle
## mle(form,start=list(h=4,log_a=2),parameters=list(log_a~poly(random)),data=rfp)

X_log_a <- model.matrix(~poly(nsize), data=rfp)

dd <- list(density=rfp$density, surv=rfp$surv,
           X_log_a=X_log_a)

obj <- MakeADFun(
  data = dd,
  parameters = list(log_a_param=c(2,1), h=4),
  DLL = "submodel")

## sdreport(obj)

opt2 <- with(obj, optim(par = par, fn = fn, gr=gr, method = "BFGS", hessian = F))

print(opt3 <- with(obj, nlminb(start = par, obj = fn, gr=gr)))


## Check
log_a <- c(X_log_a %*% c(opt3$par[1:2]))
h <- opt3$par[3]

binom <- function(size, prob, x) {log(size)+x*log(prob)+(size-x)*log((1-prob))}
binom(size = rfp$density, prob = exp(log_a)/(1 + exp(log_a)*h*rfp$density), x= rfp$surv)
## getting NAs, prob > 1
print(prob <- exp(log_a)/(1 + exp(log_a)*h*rfp$density))

log_a <- c(X_log_a %*% c(fit@coef[1:2]))
h <- fit$coef[3]
d <- dbinom(rfp$surv, size = rfp$density, prob = exp(log_a)/(1 + exp(log_a)*h*rfp$density))

## step by step check
# fun4 <- mkfun(form,start=list(h=4,log_a=2),
#               parameters=list(log_a~poly(random)),data=rfp)
# opt4 <- with(fun4, nlminb(start = unlist(start), obj = fn, gr=gr))

