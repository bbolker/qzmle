
#' @examples
#' d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
#' fit0 <- mle(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#' ## compare with bbmle
#' fit1 <- bbmle::mle2(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#'
#'fit3 <- mle(y~dnorm(mean=ymean, sd=2),start=list(ymean=mean(d$y), ysd=2),data=d)
#'fit4 <- bbmle::mle2(y~dnorm(mean=ymean, sd=2),start=list(ymean=mean(d$y),ysd=2),data=d)
#' @export
mle <- function(form, start, data, optCtrl=list(method="BFGS")) {
    ff <- mkfun(form, data)
    argList <- list(par=unlist(start), fn=ff$fn, gr=ff$gr)
    opt <- do.call(optim, c(argList,optCtrl))
    result <- list()
    result$call <- form ## need to fix to print out all input arg
    result$coefficients <- opt$par
    result$minuslogl <- opt$value
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
summary.qzmle <- function(object) {
    cat("hello\n")
}

## S3 methods for `print`
## `summary`: check the summary method for bbmle
## summary methods typically 

## `coef`, `vcov` are simpler
## (you don't really need to define the coef() method

## S3 methods in bbmle
## m <- new("mle2", call=call, call.orig=call.orig, coef=coef,
## fullcoef=unlist(fullcoef), vcov=tvcov,
## min=min, details=oout, minuslogl=minuslogl, method=method,
## optimizer=optimizer,data=as.list(data),formula=formula)
## attr(m,"df") = length(m@coef)


## twoWords, two.words, two_words, twowords, TwoWords
check_dots <- function(...) {
    if (length(list(...))>0) {
        stop("unused parameters passed to method")
    }
}
