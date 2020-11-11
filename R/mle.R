
#' @examples
#' d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
#' fit0 <- mle(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#' ## compare with bbmle
#' fit1 <- bbmle::mle2(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
#' 
mle <- function(form, start, data, optCtrl=list(method="BFGS")) {
    ff <- mkfun(form, data)
    argList<- list(par=unlist(start), fn=ff$fn, gr=ff$gr)
    opt <- do.call(optim, c(argList,optCtrl))
    return(opt)
}
