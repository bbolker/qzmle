library(Matrix)
library(TMB)
# compile("mm.cpp")
# dyn.load(dynlib("mm"))

parameters=list(log_a ~ x1 + x2 + (1|g))

dd <- data.frame(x1=rnorm(100),x2=rnorm(100),
                 g=factor(rep(1:10,each=10)))

form <- log_a ~ x1 + x2 + (1|g)

lme4::nobars(form[-2])  ## fixed
random_part <-  lme4::findbars(form[-2])
result <- lme4::mkReTrms(random_part,dd)

Z_sparse <- t(result$Zt)


