library(TMB)
library(bbmle)
compile("./submodel.cpp")
dyn.load(dynlib("./submodel"))

set.seed(101)
rfp <- transform(emdbook::ReedfrogPred,
                 nsize=as.numeric(size), random=rnorm(48))

form <- surv ~ dbinom(size = density,
                      prob = exp(log_a)/(1 + exp(log_a)*h*density))

##bbmle
mle1 <- bbmle::mle2(form,start=list(log_a=c(2,0), h=4),
                    parameters=list(log_a~poly(random)),data=rfp)


## qzmle
## mle(form,start=list(h=4,log_a=2),parameters=list(log_a~poly(random)),data=rfp)

X_log_a <- model.matrix(~poly(random), data=rfp)

dd <- list(density=rfp$density, surv=rfp$surv,
           X_log_a=X_log_a)

obj <- MakeADFun(
  data = dd,
  parameters = list(log_a_param=c(2,0), h=4),
  DLL = "submodel")

## sdreport(obj)

## opt1 uses Nelder-Mead (optim() default): it's not horribly
##  surprising that this does badly since it ignores the derivatives
##  (note that you don't get all the "outer mgc:" stuff
##  printed)
opt1 <- with(obj, optim(par = par, fn = fn, gr=gr))
opt2 <- with(obj, nlminb(start = par, obj = fn, gr=gr))

## try TMB with optim+BFGS (which *does* use derivatives)
opt3 <- with(obj, optim(par = par, fn = fn, gr=gr,
                        method="BFGS"))

## compare negative log-likelihoods
nlls <- c(TMB_NM=opt1$val,
  TMB_nlminb=opt2$objective,
  bbmle=c(-logLik(mle1)),
  TMB_BFGS=opt3$val)

## TMB with BFGS and nlminb are nearly equivalent
## (diff of 4e-10); bbmle is a little worse; TMB_NM
## is considerably worse/gets stuck
nlls-min(nlls)


## compare parameters
cbind(TMB_NM=opt1$par,TMB_nlminb=opt2$par,
      TMB_BFGS=opt3$par,bbmle=coef(mle1))

## set up calculation over the entire parameter space
ranges <- list(log_a_1=c(-1,0.5,51),
               log_a_2=c(-0.5,0.5,51),
               ## log_h = c(-8,0,9))
               log_h = c(-5,-3,11))
## vectors of values to try
## round vecs to two digits to make names pretty
vecs <- lapply(ranges,
               function(x) round(seq(x[1],x[2],length.out=x[3]),2))
names(vecs) <- names(ranges)
## set up array for results
res <- array(NA,dim=sapply(ranges,"[[",3),
             dimnames=setNames(vecs,names(ranges)))

## run (don't forget to back-transform log-h)
for (i in seq(ranges$log_a_1[3])) {
    for (j in seq(ranges$log_a_2[3])) {
        for (k in seq(ranges$log_h[3])) {
            res[i,j,k] <- obj$fn(c(vecs$log_a_1[[i]],
                                   vecs$log_a_2[[j]],
                                   exp(vecs$log_h[[k]])))
        }
    }
}

## convert to long format
df <- reshape2::melt(res)

## since panels will be values of log-h, need to find
## values closest to estimated values of log-h for each fit
get_closest <- function(x,y=unique(df$log_h)) {
    y[which.min(abs(x-y))]
}

## points corresponding to fits
pts <- as.data.frame(
    rbind(
        c(opt1$par[1:2],get_closest(log(opt1$par[3]))),
        c(opt2$par[1:2],get_closest(log(opt2$par[3]))),
        c(coef(mle1)[1:2],get_closest(log(coef(mle1)[3]))))
)
names(pts) <- c("log_a_1","log_a_2","log_h")
pts$label <- c("TMB_NM","TMB_nlminb","bbmle")


library(ggplot2); theme_set(theme_bw())
(ggplot(df, aes(log_a_1,log_a_2))
    + geom_raster(aes(fill=value-min(value,na.rm=TRUE)))
    + facet_wrap(~log_h,labeller=label_both)
    + scale_fill_viridis_c(trans="log10")
    + geom_contour(aes(z=value-min(value,na.rm=TRUE)),
                   breaks=c(1,2,10), colour="red",lty=2)
    + theme(panel.spacing=grid::unit(0,"lines"))
    + scale_x_continuous(expand=expansion(0,0))
    + scale_y_continuous(expand=expansion(0,0))
    + geom_point(data=pts, size=3,colour="red")
    + geom_label(aes(label=label),data=pts, size=3,colour="red")
)


## other checks:
## make sure that all three methods converge to the best
##   result if *started* at the best_fit result

## make sure that bbmle gives same value for best-fit
## parameter (this is where I found the problem!)

all.equal(do.call(mle1@minuslogl,as.list(unname(opt2$par))),
          opt2$objective)


## check equality of by-hand calc with TMB calc
log_a <- X_log_a %*% opt2$par[1:2]
a <- exp(log_a)
h <- opt2$par[["h"]]
prob <- a/(1+a*h*rfp$density)
R_nll <- with(rfp,
              -sum(dbinom(surv, prob=prob,
                          size=density, log=TRUE)))
all.equal(R_nll, opt2$objective)
