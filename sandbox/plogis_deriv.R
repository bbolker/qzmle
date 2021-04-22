library(Deriv)
Deriv(~x+y,"x")
Deriv(~1/(1+exp(-x))+y,"x")

f <- ~plogis(x)+y
## we COULD try to replace plogis(x) with 1/(1+exp(-x))
## worst way:
as.formula(gsub("plogis([^)]+)","abc(\\1)",deparse(f)))

drule[["plogis"]] <- alist(x=dlogis(x))
drule[["plogis"]] <- NULL
drule[["sinpi"]] <- alist(x=pi*cospi(x))
drule[["sinpi"]] <- alist(x=dlogis(x))
drule[["Plogis"]] <- alist(x=dlogis(x))
## MUST USE THE ARGUMENT AS ACTUALLY DEFINED IN THE FUNCTION
## i.e. first argument of plogis is actually called q!
drule[["plogis"]] <- alist(q=dlogis(q))
Plogis <- function(x) {}
invlogit <- function(x) { 1/(1+exp(-x))}
drule[["invlogit"]] <- alist(x=dlogis(x))
Deriv(~sinpi(x),"x")
Deriv(~Plogis(x),"x")
Deriv(~invlogit(x),"x")
Deriv(~plogis(x),"x")



                           
