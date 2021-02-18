#' @name TMB_mle_function
#' @param form A formula in expression form of "y ~ model"
#' @param parameter A list of initial values for the parameters
#' @param data A list of parameter in the formula with values in vectors
#' @param link link function and the model parameter
#' @examples
#' set.seed(123)
#' x <- runif(20, 1, 10)
#' y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
#' form <- y ~ dnorm(b0 + b1 * x, log_sigma)
#' links <- list(sigma="log")
#' TMB_mle_function(form, parameter=list(b0=0,b1=0, log_sigma=sd(y)),
#' data=list(x=x,y=y),links)

TMB_mle_function <- function(formula,parameter,data=NULL,link=NULL) {
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
  }


  ## parse distribution
  y <- formula[[2]]
  RHS <- formula[[3]]
  ddistn <- as.character(RHS[[1]])

  ## check link
  if (is.null(link)){
    link_f <- "identity"
    pname <- loglik_list[[ddistn]]$params
  } else{
    link_f <- link[[1]]
    pname <- names(link)
  }


  ## store all parameters
  params <- c()
  for (i in seq_along(parameter)) {
    params <- paste0(params, sprintf("PARAMETER(%s); ", names(parameter)[i]))
  }

  ## cpp script
  header <- "#include <TMB.hpp>

  template<class Type>
  Type objective_function<Type>::operator() () { "



  links <- paste0("Type ", invlinkfun(pname, link_f), "; ")

  nll <- sprintf("Type nll = 0.0;\nnll = -sum(%s(%s, %s, %s, true));\nreturn nll;\n}",
                 as.character(RHS[[1]]), as.character(y), toString(RHS[]), pname)

  model <- paste0(header, data_var, params, links, nll)

  write(model, file='template.cpp')
  #TMB::compile("template.cpp")
  #print("finish compiling")
}



