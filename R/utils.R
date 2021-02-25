## Link functions for TMB
## make parameter name
plinkfun <- function(pname, linkname) {
  ifelse(linkname=="identity",pname,
         paste(linkname, pname, sep="_"))
}


## names of all allowed links (except 'identity')
all_links <- c("log","logit","cloglog","sqrt","inverse","log10")
trans_parnames <- function(p) {
    regex <- sprintf("(%s)_", paste(all_links,collapse="|"))
    gsub(regex,"",p)
}




## put the rest of the pieces together ...
## then maybe filter out the identity ones
##  whatever[linkname!="identity"]

## make parameter name
invlinkfun <- function(pname, linkname) {
  switch (linkname,
    log = sprintf("%s = exp(%s)", pname, plinkfun(pname, linkname)),
    logit = sprintf("%s = invlogit(%s)", pname, plinkfun(pname, linkname))
  )
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
