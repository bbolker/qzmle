## HOW TO MAKE THIS REQUIRE PACKAGE
usethis::use_package("Deriv")



# List of log-lik function for different ditributions
loglik_list <- list(
  dpois = list(expr=expression(x * log(lambda) - lambda - lfactorial(x)),
               params=c("lambda")),
  dnorm = list(expr=expression(-log(2*4*(4*atan(1/5)-atan(1/239)))/2-log(sd) - (x-mean)^2/(2*sd^2)),
               params=c("mean","sd"))
)


#### need to convert pi as a build-in constant??? what if pi is a parameter
### take care of maybe the e constant



## data frame lat, long
y ~ dpois(exp(log_lambda), ...,
          parameters = list(log_lambda = ~ poly(lat, long, 2))
)

#' @export
mkfun <- function(formula, data) {
  if(missing(data)) {
    stop(paste("missing data...")) # if no data
    }
  RHS <- formula[[3]] # dpois(lambda = (b0 * latitude^2))
  response <- formula[[2]] #y
  ddistn <- as.character(RHS[[1]]) ## dpois /// get the name of distribution variable
  arglist <- as.list(RHS[-1]) ## $lambda = (b0 * latitude^2) ///delete function name
  arglist1 <- c(
    list(x = response),
    arglist, ##
    list(log = TRUE)
  )
  fn <- function(pars) { ## parameter
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
    d0 <- Deriv::Deriv(LL, mnames)
    ## evaluate all of the arguments to the log-likelihood
    arglist_eval <- lapply(arglist, eval, pars_and_data) ##lambda
    ## evaluate response variable and assign its value to 'x'
    arglist_eval$x <- eval(response, pars_and_data) #x = y
    ## derivative of log-lik wrt PDF parameters
    d1 <- eval(d0, arglist_eval) ## sub back to d0
    ## compute the deriv of log_lik with respect to its parameters
    parnames <- setdiff(all.vars(RHS), names(data))
    dlist <- list()
    glist <- list()
    for (p in parnames) {
      dlist[[p]] <- eval(Deriv::Deriv(arglist$lambda, p), pars_and_data)
      glist[[p]] <- -sum(d1 * dlist[[p]])
    }
    return(unlist(glist))
    ## d(loglik_pois/d(lambda))* d(lambda)/d(b0)
  }
  return(list(fn = fn, gr = gr))
}
