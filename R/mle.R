## mle function
#' Deriving MLE
#' @name mle
#' @param form A formula in expression form of "y ~ model"
#' @param start A list of initial values for p
#' @param data A list of parameter in the formula with values in vectors
#' @param optCtrl Optimization method to use
#' @examples
#' d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
#' fit0 <- mle(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#' fit3 <- mle(y~dnorm(mean=ymean, sd=2),start=list(ymean=mean(d$y), ysd=2),data=d)
#' @export

#' @importFrom stats optim
#' @importFrom numDeriv jacobian
#' @importFrom MASS ginv
mle <- function(form, start, data, optCtrl=list(method="BFGS")) {
    ff <- mkfun(form, data)
    argList <- list(par=unlist(start), fn=ff$fn, gr=ff$gr)
    opt <- do.call(stats::optim, c(argList,optCtrl))
    ## vcov is the inverse of hessian (jacobian of gradient)
    hess <- numDeriv::jacobian(ff$gr, opt$par)
    tvcov <- MASS::ginv(hess) ## solve() gives error for non-invertible matrix
    colnames(tvcov) <- names(opt$par)
    rownames(tvcov) <- names(opt$par)

    result <- list()
    result$call <- match.call()
    result$coefficients <- opt$par
    result$minuslogl <- opt$value
    result$tvcov <- tvcov
    class(result) <- "qzmle"
    return(result)
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
}


#' @export
logLik.qzmle <- function(object, ...) {
  check_dots(...)
  cat("'log Lik.'", -1*round(object$minuslogl, 2))
  cat(" (df=")
  cat(length(object$coefficients))
  cat(")")
}


#' @export
vcov.qzmle <- function(object, ...) {
  check_dots(...)
  print(object$tvcov)
}


## not working
## fit4 <- mle(y~dnorm(mean=ymean, sd=ysd),start=list(ymean=mean(d$y), ysd=sd(d$y)),data=d)


## compare with bbmle
##fit1 <- bbmle::mle2(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
##fit4 <- bbmle::mle2(y~dnorm(mean=ymean, sd=2),start=list(ymean=mean(d$y),ysd=2),data=d)

## S3 methods for `print`
## `summary`: check the summary method for bbmle
## summary methods typically

## `coef`, `vcov` are simpler
## (you don't really need to define the coef() method


## twoWords, two.words, two_words, twowords, TwoWords
check_dots <- function(...) {
    if (length(list(...))>0) {
        stop("unused parameters passed to method")
    }
}
