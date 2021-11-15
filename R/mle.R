## mle function
#' Deriving MLE
#' @name mle
#' @param form A formula in expression form of "y ~ model"
#' @param start A list of initial values for parameters
#' @param data A data frame containing any variables used in the formula for the log-likelihood. (\code{mle} does not work with variables taken from the global environment; if necessary, wrap these variables in \code{data.frame} to pass them to the function.)
#' @param fixed A list of parameter in the formula to keep fixed during optimization
#' @param parameters A list of linear submodels and random effects
#' @param control A list of parameter to pass to optimizer (See `mle_control`)
#' @param links link function for parameters (identity link as default)
#' @param method base R or TMB integration
#' @param random ???
#' @examples
#' d <- data.frame(x = 0:10, y = c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
#' fit0 <- mle(y ~ dpois(lambda = ymean), start = list(ymean = mean(d$y)), data = d)
#' ss <- list(ymean = mean(d$y), logsd = log(sd(d$y)))
#' fit3 <- mle(y ~ dnorm(mean = ymean, sd = exp(logsd)), start = ss, data = d)
#'
#' set.seed(123)
#' x <- runif(20, 1, 10)
#' y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
#' form <- y ~ dnorm(b0 + b1 * x, log_sigma)
#' fit <- mle(form,
#'   start = list(b0 = 1, b1 = 2, log_sigma = sd(y)),
#'   data = list(x = x, y = y), links = c(b0 = "identity", b1 = "identity", sigma = "log")
#' )
#' ## linear submodels
#' set.seed(101)
#' rfp <- transform(emdbook::ReedfrogPred, nsize = as.numeric(size), random = rnorm(48))
#' form <- surv ~ dbinom(size = density, prob = exp(log_a) / (1 + exp(log_a) * h * density))
#' fit4 <- mle(form, start = list(h = 4, log_a = 2),
#'             parameters = list(log_a ~ poly(random)), data = rfp)
#' @export
## bbfit4 <- bbmle::mle2(form, parameters=list(log_a~poly(random)),start=list(log_a=c(2,0), h=4),data=rfp)

#' @importFrom stats optim optimHess make.link
#' @importFrom numDeriv jacobian hessian
#' @importFrom MASS ginv
#' @importFrom lme4 findbars
mle <- function(form, start, data,
                fixed = NULL,
                links = NULL,
                parameters = NULL,
                random = NULL,
                control = mle_control(),
                method = c("R", "TMB")) {

  ## check bad link functions
  if (!is.null(links)) {
    bad_links <- which(is.na(unname(all_links[unlist(links)])))
    if (length(bad_links) > 0) {
      stop("undefined link(s): ", paste(links[bad_links], collapse = ", "))
    }

    ## Automatically translate scales
    ## FIXME: need to use mklinkfun as check
    # plinkscale <- numeric(length(start))
    # names(plinkscale) <- plinkfun(names(start), links)
    # for (i in seq_along(plinkscale)) {
    #   mm <- make.link(links[[i]])
    #   plinkscale[[i]] <- mm$linkfun(start[[i]])
    # }
  }


  ## calling TMB integration if chose to
  method <- match.arg(method)

  ## check for random effect and send to TMB
  if (!is.null(parameters)) {
    reTerms <- unlist(sapply(parameters, lme4::findbars))
    if (!is.null(reTerms)) {
      if (method != "TMB") warning("Computing random effects in TMB")
      method <- "TMB"
    }
  }


  ## make obj fun and gradients
  ff <- switch(method,
    TMB = TMB_mkfun(form, start, links, parameters, data),
    R =  mkfun(form, start, links, parameters, data),
    stop(paste("unknown method", sQuote(method)))
  )

  ## optimize using optim BFGS
  argList <- list(par = unlist(ff$start), fn = ff$fn, gr = ff$gr)
  opt <- do.call(stats::optim, c(argList, control$optControl))


  ## ------------
  ## check for fixed parameters
  if (!is.null(fixed) && !is.list(fixed)) {
    if (is.null(names(fixed)) || !is.vector(fixed)) {
      stop("'fixed' must be a named vector or named list")
    }
    fixed <- as.list(fixed)
  }
  ## FIXME: check consistency with start values

  skip_hessian <- FALSE
  if (control$hessian_method == "none") {
    skip_hessian <- TRUE
  } else {
    ##  if we don't have gradient, we'll need to use stats::optimHess() (if hessian_method=="simple")
    ##    *or* numDeriv::hessian() instead (if hessian_method=="Richardson")
    if (is.null(argList$gr)) {
      if (control$hessian_method == "simple") {
        hess <- stats::optimHess(argList)
      }
      else if (control$optControl == "Richardson") {
        hess <- numDeriv::hessian(argList$fn, opt$par, method = control$optControl)
      }
    }
    ##  when we have gradient
    ##  if the Hessian is bad AND method=="simple"
    ##    change method to "Richardson"
    else {
      hess <- numDeriv::jacobian(ff$gr, opt$par, method = control$hessian_method)
      if ((!hessian_check(hess)) && (control$hessian_method == "simple")) {
        hess <- numDeriv::jacobian(ff$gr, opt$par, method = "Richardson")
      }
    }
  }

  ## if hessian is "none"
  if (skip_hessian) {
    tvcov <- matrix(NA, length(opt$par), length(opt$par))
  } else {
    tvcov <- MASS::ginv(hess)
  }
  dimnames(tvcov) <- list(names(opt$par), names(opt$par))

  mc <- match.call()
  response <- eval(eval(mc$form)[[2]], data)
  nobs <- if (is.matrix(response)) nrow(response) else length(response)

  result <- list(
    call = mc,
    fixed = fixed,
    coefficients = opt$par,
    minuslogl = opt$value,
    formula = mc$form,
    tvcov = tvcov,
    nobs = nobs
  )
  class(result) <- "qzmle"
  return(result)
}


#' return default values and/or user-set values for details of fitting
#' @param optControl list of control parameters for \code{optim}
#' @param hessian_method method for numerically computing Hessian
#' @export
mle_control <- function(optControl = list(
                          method = "BFGS",
                          control = list()
                        ),
                        hessian_method = c("simple", "Richardson", "none")) {

  ## Don't allow other optimizer methods (yet)
  if (is.null(optControl$method) || (optControl$method != "BFGS")) {
    optControl$method <- "BFGS"
  }

  ## we want to compute our own hessian
  optControl$hessian <- FALSE

  hessian_method <- match.arg(hessian_method)
  return(named_list(optControl, hessian_method))
}



#' @export
print.qzmle <- function(x, ...) {
  check_dots(...)
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nLog-likelihood: ")
  cat(-1 * round(x$minuslogl, 2), "\n")
}


#' @export
#' @importFrom stats pnorm printCoefmat
summary.qzmle <- function(object, ...) {
  check_dots(...)
  cat("Maximum likelihood estimation\n\nCall:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  cmat <- cbind(
    Estimate = object$coefficients,
    `Std. Error` = sqrt(diag(object$tvcov))
  )
  zval <- cmat[, "Estimate"] / cmat[, "Std. Error"]
  pval <- 2 * stats::pnorm(-abs(zval))
  coefmat <- cbind(cmat, "z value" = zval, "Pr(z)" = pval)
  stats::printCoefmat(coefmat)
  cat("\n-2 log L:", 2 * object$minuslogl, "\n")
} ## maybe break it down to print.summary.qzmle?


#########################
#' @export
logLik.qzmle <- function(object, ...) {
  check_dots(...)
  ##cat("'log Lik.'", -1 * round(object$minuslogl, 2))
  ##cat(" (df=")
  ##cat(length(object$coefficients))
  ##cat(")")
  val <- -1* round(object$minuslogl, 2)
  attr(val, "df") <- length(object$coefficients)
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}


#########################
#' @export
vcov.qzmle <- function(object, ...) {
  check_dots(...)
  object$tvcov
}

## need to define stdEr generic or import/export from misctools
## #' @export
## stdEr.qzmle <- function(object, ...){
##   check_dots(...)
##   print(sqrt(diag(object$tvcov)))
## }


check_dots <- function(...) {
  if (length(list(...)) > 0) {
    stop("unused parameters passed to method")
  }
}
