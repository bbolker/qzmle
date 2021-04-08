library(TMB)
library(bbmle)
compile("./submodel.cpp")
dyn.load(dynlib("./submodel"))

set.seed(101)
rfp <- transform(emdbook::ReedfrogPred,
                 nsize=as.numeric(size), random=rnorm(48))

form <- surv ~ dbinom(size = density,
                      prob = exp(log_a)/(1 + exp(log_a)*h*density))

##bbmle
mle1 <- bbmle::mle2(form,start=list(log_a=c(2,0), h=4),
                    parameters=list(log_a~poly(random)),data=rfp)


## qzmle
## mle(form,start=list(h=4,log_a=2),
## parameters=list(log_a~poly(random)),
## data=rfp)

X_log_a <- model.matrix(~poly(random), data=rfp)

dd <- list(density=rfp$density, surv=rfp$surv,
           X_log_a=X_log_a)

obj <- MakeADFun(
  data = dd,
  parameters = list(log_a_param=c(2,0), h=4),
  DLL = "submodel")

## sdreport(obj)

opt2 <- with(obj, nlminb(start = par, obj = fn, gr=gr))

