#' @importFrom Deriv Deriv

## List of log-lik function for different distributions

loglik_list <- list(
  dpois = list(expr=expression(x * log(lambda) - lambda - lfactorial(x)),
               params=c("lambda")),
  dnorm = list(expr=expression(## -log(2*4*(4*atan(1/5)-atan(1/239)))/2 -
                   - log(2*pi)
                   - log(sd)
                   - (x-mean)^2/(2*sd^2)),
               ## this is another
               params=c("mean","sd"))
)

## data frame lat, long
y ~ dpois(exp(log_lambda), ...,
          parameters = list(log_lambda = ~ poly(lat, long, 2))
)

#' Deriving the log-lik and gradients
#' @param formula A formula in expression form of "y ~ model"
#' @param data A list of parameter in the formula with values in vectors
#' @examples
#' set.seed(101)
#' dd <- data.frame(y=rpois(100,lambda=1))
#' fun1 <- mkfun(y~dpois(exp(log_mu)), data=dd)
#' @export
mkfun <- function(formula, data) {
    ## explicit error message: otherwise won't get caught until
    ## much later
  if(missing(data)) {
    stop("missing data...") # if no data
  }
  RHS <- formula[[3]] # dnorm(mean = b0 + b1 * latitude^2, sd = 1)
  response <- formula[[2]] # always y
  ddistn <- as.character(RHS[[1]]) ## dnorm /// get the name of distribution variable
  arglist <- as.list(RHS[-1]) ## $lambda = (b0 * latitude^2), sd///delete function name
  arglist1 <- c(
    list(x = response), ##assign x to y)
    arglist,
    list(log = TRUE)
  )
  fn <- function(pars) { ## parameter
    pars_and_data <- c(as.list(pars), data) ## list of b0,b1,y,lattitude
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

    parnames <- setdiff(all.vars(RHS), names(data))
    glist <- list()
    ## a matrix with appropriately named columns corresponding to parameters
    ## we  know what the structure of the returned gradient vector (which
    ## constitutes the gradients for all observations squashed together,
    ## i.e.  g_11, g_12,... g_21,g_22,... g_ij
    ##  where i indicates the parameter and j indicates the observation
    ## uses the fact that R stores matrix in column-major order
    d1_mat <- matrix(d1, ncol=length(mnames), dimnames=list(NULL, mnames))
    for (m in mnames){
        d2 <- d1_mat[,m]
        if (is.numeric(arglist[[m]])) {
            ## constant!
            glist[[m]] <- 0
        } else {
            for (p in parnames){
                dlist <- list()
                dlist[[m]][[p]] <- eval(Deriv::Deriv(arglist[[m]],p), pars_and_data)
                glist[[m]][[p]] <- -sum(d2*dlist[[m]][[p]])
            } ## p in parnames
        } ## arg is not constant
    } ## m in mnames
    return(unlist(glist))

    ##sd - d(loglik_norm)/d(sd)
    ##b0 - d(loglik_norm)/d(norm) * d(mean)/d(b0)
    ##b1 - d(loglik_norm)/d(norm) * d(mean)/d(b1)

    ## d(loglik_pois/d(lambda))* d(lambda)/d(b0)
  }
  return(list(fn = fn, gr = gr))
}
