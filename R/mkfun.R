## List of log-lik function for different distributions
loglik_list <- list(
  dbinom = list(expr=expression(log(size)+x*log(prob)+(size-x)*log((1-prob))),
                params=c("prob")),
  dpois = list(expr=expression(x * log(lambda) - lambda - lfactorial(x)),
               params=c("lambda")),
  dnorm = list(expr=expression(
    - log((2*pi)^0.5)
    - log(sd)
    - (x-mean)^2/(2*sd^2)),
    params=c("mean","sd"))
)

#' Allow users to add their own log likelihood functions
#' @name add_logl
#' @param funct objective function (must be defined in environment)
#' @param logl log likelihood function as expression
#' @param params a vector of parameters for the objective function
#' @export
add_logl <- function(funct, logl, params){
  if(try(check_fun(funct))){
    funct_name <- as.character(substitute(funct))
    new_funct <- list(expr=logl, params=params)
    new_list <- list(new_funct)
    names(new_list) <- funct_name
    loglik_list <- c(loglik_list, new_list)
  }
}

#' Deriving the log-lik and gradients
#' @name mkfun
#' @param formula A formula in expression form of "y ~ model"
#' @param data A list of parameter in the formula with values in vectors
#' @param links Link function for each parameters
#' @param parameters A list of linear submodels
#' @examples
#' set.seed(101)
#' dd <- data.frame(y=rpois(100,lambda=1))
#' fun1 <- mkfun(y~dpois(exp(lambda)), data=dd)
#' fun2 <- mkfun(y~dnorm(mean = b0 + b1 * latitude^2, sd = 1), data=dd)
#' form <- surv ~ dbinom(size = density, prob = 1/(1 + exp(log_a)*h*density))
#' fun3 <- mkfun(form,parameters=list(log_a~poly(size)),data=emdbook::ReedfrogPred)
#' @export

## trying to replicate:
## library(emdbook)
## bbmle::mle2(surv ~ dbinom(size=density,prob=1/(1+exp(log_a)*h*density)),
##     parameters=list(log_a~size),
##     start=list(log_a=0,h=1),
##     data=ReedfrogPred)



#' @importFrom Deriv Deriv
mkfun <- function(formula,
                  links=NULL,
                  parameters=NULL,
                  data) {
  if(missing(data)) {
    stop("missing data...") # if no data
  }
  RHS <- formula[[3]]
  response <- formula[[2]]
  ddistn <- as.character(RHS[[1]])

  ## Check distribution
  ## suggest to add user's own likelihood function
  if(try(check_fun(ddistn))){
    if (!ddistn %in% names(loglik_list)) {
      stop("Can't evaluate the likelihood for ", sQuote(ddistn),
           paste("\n Use add_logl() to add the log likelihood function"))
    }
  }

  ## submodels
  if(!missing(parameters)){
    ##FIXME: assumed only one submodel for now
    submodel_vars <- as.character(sapply(parameters,"[[",2))
    dvars <- parameter_parse(parameters, data)
  } else{
      submodel_vars = NULL
    }

  ## assign distribution parameters
  mnames <- loglik_list[[ddistn]]$params
  arglist <- as.list(RHS[-1]) ## $lambda = (b0 * latitude^2), sd///delete function name

  ## in case user puts other parnames ie. mu instead of mean
  ## names(arglist) <- mnames

  arglist1 <- c(
    list(x = response), ##assign x to y)
    arglist,
    list(log = TRUE)
  )

  ## do we want likelihood respect to orig or link
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

    ## eventually we need to calculate partial derivatives of the log-likelihood
    ## with respect to all of its parameters
    LL <- loglik_list[[ddistn]]$expr
    d0 <- Deriv::Deriv(LL, mnames) ## d(dist)/d(mnames)
    arglist_eval <- lapply(arglist, eval, pars_and_data) ##mean, sd
    arglist_eval$x <- eval(response, pars_and_data) ##evaluate response variable and assign its value to 'x'
    d1 <- eval(d0, arglist_eval) ## sub d0 - compute the deriv of log_lik wrt to its parameters
    ## d1 = D(dbinom/prob)

    ## parameters of model parameter
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
              if(p %in% all.vars(arglist[[m]])) {
                ## links
                # tlink <- links[[p]]
                # mm <- make.link(tlink)

                dlist <- list()
                ## d(mean)/d(b0); d(prob)/d(log_a)
                dlist[[m]][[p]] <- eval(Deriv::Deriv(arglist[[m]], p), pars_and_data)

                ## check if parameter has submodel
                if(p %in% submodel_vars){
                  vars <- names(dvars)
                  for (i in seq_along(vars)){
                    ## d(log_a)/d(log_a.intercept)
                    d3 <- dvars[[i]]*(i-1)
                    glist[[m]][[vars[i]]] <- -sum(d3*d2*dlist[[m]][[p]])
                  }
                } else{
                  # deriv rule on links - d(b0)/d(log_b0)
                  ##dlist[[m]][[p]] <- 1/mm$mu.eta(mm$linkinv(dlist[[m]][[p]]))
                  glist[[m]][[p]] <- -sum(d2*dlist[[m]][[p]])
                }
              }
            } ## p in parnames
        } ## arg is not constant
    } ## m in mnames
    return(unlist(glist))

    ##sd - d(loglik_norm)/d(sd) * d(sd)/d(log_sd)
    ##b0 - d(loglik_norm)/d(norm) * d(mean)/d(b0) * d(b0)/d(log(b0))
    ##b1 - d(loglik_norm)/d(norm) * d(mean)/d(b1)

    ## d(loglik_pois/d(lambda))* d(lambda)/d(b0)
  }
  return(list(fn = fn, gr = gr))
}

## setting up submodels
parameter_parse <- function(parameters, data){
  ## assume only one submodel for now
  formula <- parameters[[1]]
  var <- formula[[2]]
  RHS <- formula[-2] ## "~ size"

  ## set up model
  X <- model.matrix(RHS, data=data)

  ## convert X into a list
  val_list <- split(X, rep(1:ncol(X), each = nrow(X)))
  names(val_list) <- paste(var, colnames(X), sep=".")
  return(val_list)
}
