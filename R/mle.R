
#' @examples
#' d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
#' fit0 <- mle(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#' ## compare with bbmle
#' fit1 <- bbmle::mle2(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#'
#'fit3 <- mle(y~dnorm(mean=ymean, sd=2),start=list(ymean=mean(d$y), ysd=2),data=d)
#'fit4 <- bbmle::mle2(y~dnorm(mean=ymean, sd=2),start=list(ymean=mean(d$y),ysd=2),data=d)
mle <- function(form, start, data, optCtrl=list(method="BFGS")) {
    ff <- mkfun(form, data)
    argList<- list(par=unlist(start), fn=ff$fn, gr=ff$gr)
    opt <- do.call(optim, c(argList,optCtrl))
    ## trying to make it looks like bbmle output
    opt$Call <- form ## need to fix to print out all input arg
    opt$Coefficients <- opt$par
    opt$`Log-likelihood` <- opt$value
    output <- tail(opt, 3)
    return(output)
}

## S3 methods for `print`, `coef`, `vcov`

## S3 methods in bbmle
## m <- new("mle2", call=call, call.orig=call.orig, coef=coef,
## fullcoef=unlist(fullcoef), vcov=tvcov,
## min=min, details=oout, minuslogl=minuslogl, method=method,
## optimizer=optimizer,data=as.list(data),formula=formula)
## attr(m,"df") = length(m@coef)




