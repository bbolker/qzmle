# qzmle: enhanced maximum likelihood estimation in R/TMB

The goals of this project/package are to create a new version of the `bbmle` package with the following enhancements (roughly in order of importance):

- choice between an R backend and a TMB backend; the latter will allow random effects in the model
- R backend also does automatic differentiation (via the `Deriv` package) when using the formula interface
- TMB backend will enable HMC sampling of the model via the `tmbstan` package
- general cleanup/reorganization
- for fixing parameters, substitute `map` (as in the TMB package) for `fixed`
- ability to add priors/regularizing factors

This package started as an undergraduate thesis by Queenie Zheng (hence the name) at McMaster University.

More details will eventually be available in a vignette.
