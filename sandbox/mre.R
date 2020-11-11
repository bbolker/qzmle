devtools::load_all()
library(Deriv)

form1 <- y ~ dpois(lambda = b0 + b1 * latitude^2)
form2 <- y ~ dnorm(mean = b0 + b1 * latitude^2, sd = 1)

set.seed(101)
dd <- data.frame(
  y = rpois(100, lambda = 2),
  latitude = rnorm(100))

pars_and_data <- as.list(c(as.list(c(b0=2, b1=3)), dd))

pois <- loglik_list$dpois$expr
pois_names <- loglik_list$dpois$params
arglist1 <- c(list(y = form1[[2]]), as.list(form1[[3]][-1]))
arglist_eval1 <- lapply(arglist1, eval, pars_and_data)
d0 <- Deriv::Deriv(pois, pois_names)
d1 <- eval(d0, arglist_eval)


normal <- loglik_list$dnorm$expr
norm_names <- loglik_list$dnorm$params
arglist2 <- c(list(y = form2[[2]]), as.list(form2[[3]][-1]))
arglist_eval2 <- lapply(arglist2, eval, pars_and_data)
d2 <- Deriv::Deriv(normal, norm_names)
d3 <- eval(d2, arglist_eval2)

