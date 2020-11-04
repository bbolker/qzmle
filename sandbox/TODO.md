## to do (in no particular order)

- make a `qzmle` repository that contains the skeleton of a `qzmle` R package that will *eventually* be a cleaner version (replacement? refactoring tool?) for `bbmle`
    - `devtools`, `usethis`
	- copy `mkfun` from `sandbox/tmp.R`: **think about a good name for it**
	- think about conventions: how are you going name functions? what code style are you going to use? https://style.tidyverse.org/ (`styler` package?)
	
- extend `mkfun` to deal with more parameters
	- what if PDFs depend on multiple parameters?
    - what if PDF parameters depend on multiple model parameters?
    - for a PDF with multiple parameters:
	     - write out the chain-rule process by hand and see how it generalizes to >1 parameter
		 - for 1-parameter: dL/d(a,b,c) = dL/d(lambda) * d(lambda)/d(a,b,c) = ({\partial{lambda}/\partial{a}, \partial{lambda}/\partial{b}, \partial{lambda}/\partial{c}})
		 - {mean, sd}  {a, b, c}:  L <- {mean, sd}, mean -> {a,b} sd -> {c}

Which parts of a formula are parameters? Which parts are data?

we have p  PDF parameters (e.g. for Poisson, p=1 {lambda})
and for each parameter we have p_i model parameters (e.g. b0+b1*latitude^2, p_1=2)

we have N observations, so dL/dlambda = vector of length N
d lambda/d(b0) = vector of length N
d lambda/d(b1) = vector of length N
