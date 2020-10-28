
devtools::load_all()
library(Deriv)
ee <- loglik_list$dnorm$expr
Deriv(ee,"mean")
Deriv(ee,"sd")
mkfun(ee)
