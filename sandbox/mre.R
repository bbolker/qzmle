
devtools::load_all()
library(tidyverse)
tilibrary(Deriv)
ee <- loglik_list$dnorm$expr
Deriv(ee,"mean")
Deriv(ee,"sd")
mkfun(ee)
