#' @name TMB_mle_function
#' @param form A formula in expression form of "y ~ model"
#' @param parameter A list of initial values for the parameters
#' @param data A list of parameter in the formula with values in vectors
#' @examples
#' set.seed(123)
#' x <- runif(20, 1, 10)
#' y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
#' form <- y ~ dnorm(b0 + b1 * x, sd)
#' TMB_mle_function(form, parameter=list(b0=0,b1=0, sd=sd(y)), data = list(x=x,y=y))

TMB_mle_function <- function(formula,
                             parameter,
                             data=NULL) {
  if (!is.null(data)) {
    data_var <- character(length(data))

    ## check data
    for (i in seq_along(data)) {
        if (is.character(data[[i]])) {
            stop("Cannot process string data")
        }
        ## store data
        data_var[i] <- sprintf("DATA_VECTOR(%s);", names(data)[i])
    }
    data_var <- paste(data_var, collapse="\n")
  }

  ## parse distribution
  y <- formula[[2]]
  RHS <- formula[[3]]
  ddistn <- as.character(RHS[[1]])

  ## store all parameters
  params <- c()
  for (i in seq_along(parameter)) {
    params <- paste0(params, sprintf("PARAMETER(%s); ", names(parameter)[i]))
  }

  ##FIX.ME: link function


  ## cpp script
  header <- "#include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() () { "

  nll <- sprintf("Type nll = 0.0;\nnll = -sum(%s(%s, %s, true));\nreturn nll;\n}",
                 as.character(RHS[[1]]), as.character(y), toString(RHS[-1])) ## not stable

  model <- paste0(header, data_var, params, nll)

  write(model, file='template.cpp')
  TMB::compile("template.cpp")
  print("finish compiling")
}



