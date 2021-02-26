## list of link to inverse link functions
all_links <- c("logit"="invlogit(%s)",
               "probit"="pnorm(%s)",
               "cauchit"=NA,
               "cloglog"=NA,
               "identity"="%s",
               "log"="exp(%s)",
               "sqrt"="%s**2",
               "1/mu^2"="1/sqrt(%s)",
               "inverse"="(1/%s)")

## make parameter name
plinkfun <- function(pname, linkname) {
  ifelse(linkname=="identity",pname,
         paste(linkname, pname, sep="_"))
}

## List of log-lik function for different distributions
loglik_list <- list(
  dpois = list(expr=expression(x * log(lambda) - lambda - lfactorial(x)),
               params=c("lambda")),
  dnorm = list(expr=expression(## -log(2*4*(4*atan(1/5)-atan(1/239)))/2 -
                   - log((2*pi)^0.5)
                   - log(sd)
                   - (x-mean)^2/(2*sd^2)),
               params=c("mean","sd"))
)

## copied from lme4
named_list <- function (...)
{
    L <- list(...)
    snm <- sapply(substitute(list(...)), deparse)[-1]
    if (is.null(nm <- names(L)))
        nm <- snm
    if (any(nonames <- nm == ""))
        nm[nonames] <- snm[nonames]
    setNames(L, nm)
}

setNames <- function (object = nm, nm)
{
  names(object) <- nm
  object
}

## check if hessian is positive definite
hessian_check <- function(x, tol=1e-08) {
    eigenvalues <- eigen(x, only.values = TRUE)$values
    n <- nrow(x)
    for (i in 1:n) {
        if (abs(eigenvalues[i]) < tol) {
            eigenvalues[i] <- 0
        }
    }
    if (any(eigenvalues <= 0)) {
        return(FALSE)
    }
    return(TRUE)
}
