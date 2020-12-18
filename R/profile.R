#' @name profilefun
#' @param fitted qzmle object
#' @param which list of parameters to profile (default: all parameters)
#' @param maxsteps number of steps to take looking for zmax
#' @param alpha max alpha level
#' @param zmax max log-likelihood difference to search to
#' @param del stepsize

profilefun <- function (fitted, which = 1:p, maxsteps = 100,
                        alpha=0.1, zmax=sqrt(qchisq(1 - alpha/2, p)),
                        del = zmax/5, ...) {
  Pnames <- names(B0 <- fitted$coefficients)
  p <- length(Pnames)
  std.err <- sqrt(diag(fitted$tvcov))
  fix0 <- fitted$fixed ## fixed parameters
  indx_fixed <- seq_along(fix0)[!is.na(fix0)] ## index of non-fixed parameters

  onestep <- function(step){
    ## sgn is c(-1,1)
    ## i-th parameter
    bi <- B0[i] + sgn * step * del * std.err[i] ## parameter going forward/backward in dx*previous step*std dev scale

  }
}
