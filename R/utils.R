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


## used in mle.R
## make param name with link function prefix
plinkfun <- function(pname, linkname) {
  ifelse(linkname=="identity",pname,
         paste(linkname, pname, sep="_"))
}


## used in TMB.R
## get param names from "linkfun_pname"
trans_parnames <- function(p) {
  regex <- sprintf("(%s)_", paste(names(all_links),collapse="|"))
  gsub(regex,"",p)
}

## submodels
parameter_parse <- function(formula, data){
  X <- model.matrix(formula, data=data)
  return(X)
}

re_parameter_parse <- function(formula, data){
  rt <- lme4::mkReTrms(formula, fr=data)
  Z <- Matrix::t(rt$Zt)
  return(Z)
}

## Check if objective function exists
## examples
## check_fun("dbinom")
## try(check_fun("djunk"))
## djunk <- function(y) {}
## try(check_fun("djunk"))
## rm(djunk)
check_fun <- function(f) {
    if (!exists(f)) stop("function: ",
                          paste0(sQuote(f), " doesn't exist"))
    ff <- formals(get(f))
    if (names(ff)[1]!="x") stop("first argument should be 'x'")
    if (!"log" %in% names(ff)) stop("function should have a 'log' argument")
    return(TRUE)
    }


## checks whether name starts with either a specific link
## (or if NULL) starts with "some link"_
## returns either "response" or "link"
## @example
## detect_scale(c(log_a=1)) # "link"
## detect_scale(c(some_x=1)) # "response"
## detect_scale(c(log_a=1, some_x=1)) # "link" "response"
detect_scale <- function(x, link=NULL) {
  links <- names(x)
  regex <- sprintf("(%s)_", paste(names(all_links),collapse="|"))
  ifelse(length(grep(regex, links))>0,
         return("link"),
         return("response"))

  ## output a vector of link functions
  # s <- character(length(links))
  # for (i in seq_along(links)){
  #   ifelse(length(grep(regex, links[i]))>0,
  #          s[i]<- "link",
  #          s[i]<- "response")
  # }
  # return(s)
}



## 1. check names, warn if necessary
## (i.e. warn if (direction=="invlink" && name doesn't start with linkname)
## OR (direction=="link" && name starts with linkname)
## 2. transform x (with $invlink or $linkfun)
## 3. modify names(x) appropriately (add or subtract <linkname>_)
## examples:
## invlink <- mklinkfun("log")
## x <- c(log_a=1)
## invlink(x)

mklinkfun <- function(linkname, direction=c("linkfun","linkinv")) {
  direction <- match.arg(direction)
  m <- make.link(linkname)
  f <-function(x) {
    scale <- detect_scale(x)

    if (direction=="linkfun"){
      if (scale=="link"){
        warning("applying link functions: ",
                paste(dQuote(linkname), "to possible link-scaled parameter:",
                      dQuote(names(x)), sep = " "))}
      x_name <- paste0(linkname,"_", names(x))
    }

    if (direction=="linkinv"){
      if (scale=="response"){
        warning("applying invlink functions: ",
                paste(dQuote(linkname), "to possible not link-scaled parameter:",
                      dQuote(names(x)), sep = " "))}
      x_name <- trans_parnames(names(x))
    }

    result <- m[[direction]](x)
    names(result) <- x_name
    return(result)
  }
  return(f)
}


## named variable val with variable name in a list (copied from lme4)
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


## set name to object (copied from stats)
setNames <- function (object = name, name)
{
  names(object) <- name
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

## formula utilities
LHS_to_char <- function(f) as.character(f[[2]])
RHS <- function(f) {
    stopifnot(identical(f[[1]],quote(`~`))) ## should be a formula
    return(f[[length(f)]])
}
onesided_formula <- function(f) {
    stopifnot(identical(f[[1]],quote(`~`))) ## should be a formula
    if (length(f)==2) return(f) ## already one-sided
    return(f[-2])
}
