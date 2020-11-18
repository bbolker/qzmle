library(bbmle)

output <- list(Call=1, Coefficient=2, value=3)
print(output)

d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
fit0 <- qzmle::mle(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
print(fit0)
class(fit0) <- "qzmle"
print.qzmle <- function(x) {
    cat("ha ha I won't actually show you anything\n")
}
print(fit0)
## if you want to see what's inside:
str(fit0)
unclass(fit0)
names(fit0)


## SUGGESTIONS:
##  $Call -> $call
##  $Coefficients -> coefficients
##  
library(bbmle)
fit0B <- bbmle::mle2(y~dpois(lambda=ymean),start=list(ymean=mean(d$y)),data=d)
print(fit0B) ## calling the show() method

getMethod("show", sig="mle2")
getMethod("summary", sig="mle2")
getMethod("show", sig="summary.mle2")

## define summary.qzmle -> returns an object of class 'summary.qzmle'
## define print.summary.qzmle -> prints that thing

## vcov is simpler
## coef is simpler

## experimenting with summary, class, etc.
m1 <- lm(y~x, data=d)
coef(m1) ## uses stats:::coef.default
print(s <- summary(m1)) ## uses stats:::summary.lm *and* stats:::print.summary.lm
coef(s)
## what to use to print a coefficient table prettily
printCoefmat(coef(s))
options(show.signif.stars=FALSE)
getMethod("summary",sig="mle2")
getMethod("show",sig="summary.mle2")
