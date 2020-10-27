library(Deriv)


# List of log-lik function for different ditributions
loglik_list <- list(
  dpois = list(expr=expression(x * log(lambda) - lambda - lfactorial(x)),
               params=c("lambda")),
  dnorm = list(expression(-log(2*pi*sd) - (x-mean)^2/(2*sd^2)),
               params=c("mean","sd"))
)


# for testing
form <- y ~ dpois(lambda = b1 * latitude^2)
set.seed(101)
dd <- data.frame(
  y = rpois(100, lambda = 2),
  latitude = rnorm(100)
)



## form -> objective function (i.e. a function that computes -sum(dpois(y, lambda=b0*latitude^2))
##  where the variables in the formula are *evaluated* in an environment that
##  includes the current parameter
##  values and information stored in a 'data' variable


## data frame lat, long
y ~ dpois(exp(log_lambda), ...,
          parameters = list(log_lambda = ~ poly(lat, long, 2))
)

mkfun <- function(formula, data) {
  RHS <- formula[[3]]
  response <- formula[[2]]
  ddistn <- as.character(RHS[[1]]) ## get the name of distribution variable
  arglist <- as.list(RHS[-1]) ## delete function name
  arglist1 <- c(
    list(x = response),
    arglist, ##
    list(log = TRUE)
  )
  fn <- function(pars) {
    pars_and_data <- c(as.list(pars), data)
    r <- with(
      pars_and_data,
      -sum(do.call(ddistn, arglist1))
    )
    return(r)
  }
  gr <- function(pars) {
    pars_and_data <- c(as.list(pars), data)
    if (!ddistn %in% names(loglik_list)) {
      stop("I can't evaluate the derivative for ", sQuote(ddistn))
    }
    ## eventually we need to calculate partial derivatives of the log-likelihood
    ## with respect to all of its parameters
    LL <- loglik_list[[ddistn]]$expr
    mnames <- loglik_list[[ddistn]]$params
    ## ???
    ## setdiff(all.vars(LL), "x")  ## response var should be the only non-parameter
    d0 <- Deriv(LL, mnames)
    ## evaluate all of the arguments to the log-likelihood
    arglist_eval <- lapply(arglist, eval, pars_and_data)
    ## evaluate response variable and assign its value to 'x'
    arglist_eval$x <- eval(response, pars_and_data)
    ## derivative of log-lik wrt PDF parameters
    d1 <- eval(d0, arglist_eval)
    ## compute the deriv of log_lik with respect to its parameters
    parnames <- setdiff(all.vars(RHS), names(data))
    dlist <- list()
    glist <- list()
    for (p in parnames) {
      dlist[[p]] <- eval(Deriv(arglist$lambda, p), pars_and_data)
      glist[[p]] <- -sum(d1 * dlist[[p]])
    }
    return(unlist(glist))
    ## d(loglik_pois/d(lambda))* d(lambda)/d(b0)
  }
  return(list(fn = fn, gr = gr))
}
