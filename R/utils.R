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

## get pnames from "linkfun_pname"
trans_parnames <- function(p) {
  regex <- sprintf("(%s)_", paste(names(all_links),collapse="|"))
  gsub(regex,"",p)
}


## Checks objective function
#' @examples
#' check_fun(dbinom)
#' try(check_fun(djunk))
#' djunk <- function(y) {}
#' try(check_fun(djunk))
#' rm(djunk)

check_fun <- function(f) {
    if (!exists(deparse(substitute(f)))) stop("function: ",
                                              paste0(sQuote(f), " doesn't exist"))
    ff <- formals(f)
    if (names(ff)[1]!="x") stop("first argument should be 'x'")
    if (!"log" %in% names(ff)) stop("function should have a 'log' argument")
    return(TRUE)
}


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
