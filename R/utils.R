## copied from lme4
#'
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
