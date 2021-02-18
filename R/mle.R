## mle function
#' Deriving MLE
#' @name mle
#' @param form A formula in expression form of "y ~ model"
#' @param start A list of initial values for p
#' @param data A list of parameter in the formula with values in vectors
#' @param fixed A list of parameter in the formula to keep fixed during optimization
#' @param control A list of parameter to pass to optimizer (See `mle_control`)
#' @param link link function for parameters (identity link as default)
#' @examples
#' d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
#' fit0 <- mle(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#' ss <- list(ymean=mean(d$y), logsd=log(sd(d$y)))
#' fit3 <- mle(y~dnorm(mean=ymean, sd=exp(logsd)),start=ss,data=d)
#' @export

#' @importFrom stats optim optimHess
#' @importFrom stats make.link
#' @importFrom numDeriv jacobian hessian
#' @importFrom MASS ginv
mle <- function(form, start, data, fixed=NULL, control=mle_control(),
                link='identity') {
    ff <- mkfun(form, data, link)


    ## check for fixed parameters
    if (!is.null(fixed) && !is.list(fixed)) {
      if (is.null(names(fixed)) || !is.vector(fixed)){
        stop("'fixed' must be a named vector or named list")}
      fixed <- as.list(fixed)
    }
    ## FIXME: check consistency with start values

    argList <- list(par=unlist(start), fn=ff$fn,
                    gr=ff$gr)
    opt <- do.call(stats::optim, c(argList,control$optControl))

    skip_hessian <- FALSE
    if (control$hessian_method=="none"){
      skip_hessian <- TRUE
    } else {
        ##  if we don't have gradient, we'll need to use stats::optimHess() (if hessian_method=="simple")
        ##    *or* numDeriv::hessian() instead (if hessian_method=="Richardson")
        if (is.null(argList$gr)) {
        if(control$hessian_method=="simple"){
          hess <- stats::optimHess(argList)
        }
        else if (control$optControl=="Richardson"){
          hess <- numDeriv::hessian(argList$fn, opt$par, method=control$optControl)
        }
      }
      ##  when we have gradient
      ##  if the Hessian is bad AND method=="simple"
      ##    change method to "Richardson"
      else{
        hess <- numDeriv::jacobian(ff$gr, opt$par, method=control$hessian_method)
        if ((!hessian_check(hess)) && (control$hessian_method=="simple")) {
        hess <- numDeriv::jacobian(ff$gr, opt$par, method="Richardson")
        }
      }
    }

    ## if hessian is "none"
    if (skip_hessian) {
      tvcov <- matrix(NA, length(opt$par), length(opt$par))
    } else{
      tvcov <- MASS::ginv(hess)
    }
    dimnames(tvcov) <- list(names(opt$par), names(opt$par))

    result <- list(
        call = match.call(),
        fixed = fixed,
        coefficients = opt$par,
        minuslogl = opt$value,
        tvcov = tvcov
    )
    class(result) <- "qzmle"
    return(result)
}


#' return default values and/or user-set values for details of fitting
#' @param optControl list of control parameters for \code{optim}
#' @param hessian_method method for numerically computing Hessian
#' @export
mle_control <- function(optControl=list(method="BFGS",
                                        control=list(),
                                        ...),
                        hessian_method=c("simple","Richardson","none")) {

  ## Don't allow other optimizer methods (yet)
  if (is.null(optControl$method) || (optControl$method!='BFGS')){
      optControl$method <- 'BFGS'
  }

  ## we want to compute our own hessian
  optControl$hessian <- FALSE

  hessian_method <- match.arg(hessian_method)
  return(named_list(optControl, hessian_method))
}



#' @export
print.qzmle <- function (x, ...) {
    check_dots(...)
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nLog-likelihood: ")
    cat(-1*round(x$minuslogl, 2), "\n")
}


#' @export
#' @importFrom stats pnorm printCoefmat
summary.qzmle <- function(object, ...) {
  check_dots(...)
  cat("Maximum likelihood estimation\n\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  cmat <- cbind(Estimate=object$coefficients,
                 `Std. Error`=sqrt(diag(object$tvcov)))
  zval <- cmat[, "Estimate"]/cmat[, "Std. Error"]
  pval <- 2*stats::pnorm(-abs(zval))
  coefmat <- cbind(cmat,"z value"=zval,"Pr(z)"=pval)
  stats::printCoefmat(coefmat)
  cat("\n-2 log L:", 2*object$minuslogl,"\n")
} ## maybe break it down to print.summary.qzmle?

#########################
#' @export
logLik.qzmle <- function(object, ...) {
  check_dots(...)
  cat("'log Lik.'", -1*round(object$minuslogl, 2))
  cat(" (df=")
  cat(length(object$coefficients))
  cat(")")
}


#########################
#' @export
vcov.qzmle <- function(object, ...) {
  check_dots(...)
  print(object$tvcov)
}

## need to define stdEr generic or import/export from misctools
## #' @export
## stdEr.qzmle <- function(object, ...){
##   check_dots(...)
##   print(sqrt(diag(object$tvcov)))
## }





check_dots <- function(...) {
    if (length(list(...))>0) {
        stop("unused parameters passed to method")
    }
}
