library(TMB)
compile("reg.cpp")
dyn.load(dynlib("reg"))

set.seed(123)
x <- runif(20, 1, 10)
y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
plot(x, y)

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
