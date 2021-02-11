library(TMB)
compile("reg.cpp")
dyn.load(dynlib("reg"))

set.seed(123)
x <- runif(20, 1, 10)
y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
##plot(x, y)

obj <- MakeADFun(
  data = list(x = x, y = y),
  parameters = list(b0 = 0, b1 = 0, log_sigma = 0),
  DLL = "reg")


rep <- sdreport(obj)
rep ## check estimates

lm(y~x) ## check if its the same as lm

## check if its the same as optim and nlminb
opt <- optim(par = obj$par, fn = obj$fn, gr=obj$gr)
opt2 <- nlminb(start = obj$par, obj = obj$fn, gr = obj$gr)


compile("reg_vec.cpp")
dyn.load(dynlib("reg_vec"))


obj2 <- MakeADFun(
  data = list(x = x, y = y),
  parameters = list(b0 = 0, b1 = 0, log_sigma = 0),
  DLL = "reg_vec")


with(obj2, nlminb(start = par, obj = fn, gr=gr))

dnbinom(x=0:5, mu=1:6, size=rep(0.2,6), log=TRUE)


mle2(y ~ dnorm(b0+b1*x, exp(log_sigma)),
     start=list(b0=0,b1=0,log_sigma=0),
     data=list(x=rnorm(5),y=rnorm(5)))
## -> create reg_vec.cpp
## compile it
## optimize it
