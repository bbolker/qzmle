## MAKE THIS REQUIRE PACKAGE
## usethis::use_package("Deriv")

# List of log-lik function for different ditributions
loglik_list <- list(
  dpois = list(expr=expression(x * log(lambda) - lambda - lfactorial(x)),
               params=c("lambda")),
  dnorm = list(expr=expression(-log(2*4*(4*atan(1/5)-atan(1/239)))/2-log(sd) - (x-mean)^2/(2*sd^2)),
               params=c("mean","sd"))
)



## data frame lat, long
y ~ dpois(exp(log_lambda), ...,
          parameters = list(log_lambda = ~ poly(lat, long, 2))
)

#' @export
mkfun <- function(formula, data) {
  if(missing(data)) {
    stop(paste("missing data...")) # if no data
    }
  RHS <- formula[[3]] # dnorm(mean = b0 + b1 * latitude^2, sd = 1)
  response <- formula[[2]] # always y
  ddistn <- as.character(RHS[[1]]) ## dnorm /// get the name of distribution variable
  arglist <- as.list(RHS[-1]) ## $lambda = (b0 * latitude^2), sd///delete function name
  arglist1 <- c(
    list(x = response), ##assign x to y
    arglist,
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
    ## setdiff(all.vars(LL), "x")  ## response var should be the only non-parameter
    d0 <- Deriv::Deriv(LL, mnames) ## evaluate all of the arguments to the log-likelihood
    arglist_eval <- lapply(arglist, eval, pars_and_data) ##mean, sd
    arglist_eval$x <- eval(response, pars_and_data) ##evaluate response variable and assign its value to 'x'
    d1 <- eval(d0, arglist_eval) ## sub d0 - compute the deriv of log_lik wrt to its parameters

    if (length(setdiff(names(arglist_eval), "x")) < 2) { ##one parameter distribution
    parnames <- setdiff(all.vars(RHS), names(data))
    dlist <- list()
    glist <- list()
    for (p in parnames) {
      dlist[[p]] <- eval(Deriv::Deriv(arglist$lambda, p), pars_and_data)
      glist[[p]] <- -sum(d1 * dlist[[p]])
    }
    }
    else{
      ##having more than one parameters
      for (m in mnames){
          ## mlist[[m]] <- grepl(...)
      }

    }
    return(unlist(glist))

    ##sd - d(loglik_norm)/d(sd)
    ##b0 - d(loglik_norm)/d(norm) * d(mean)/d(b0)
    ##b1 - d(loglik_norm)/d(norm) * d(mean)/d(b1)

    ## d(loglik_pois/d(lambda))* d(lambda)/d(b0)
  }
  return(list(fn = fn, gr = gr))
}
