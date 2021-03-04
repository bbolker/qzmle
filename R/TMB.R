#' @name TMB_template
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
#' TMB_template(form, parameter=list(b0=1,b1=2, log_sigma=sd(y)),
#' data=list(x=x,y=y),link=c(b0="identity", b1="identity", sigma="log"))
#' ff <- TMB_mkfun(form, parameter=list(b0=1,b1=2, log_sigma=sd(y)),
#' data=list(x=x,y=y),link=c(b0="identity", b1="identity", sigma="log"))

TMB_template <- function(formula,parameter,data=NULL,link=NULL) {
  if (!is.null(data)) {
    data_var <- character(length(data))

    ## check data
    for (i in seq_along(data)) {
        if (is.character(data[[i]])) {
            stop("Cannot process string data")
        }
        ## store data
        data_var[i] <- sprintf("DATA_VECTOR(%s); \n", names(data)[i])
    }
  }

  ## parse distribution
  y <- formula[[2]]
  RHS <- formula[[3]]
  ddistn <- as.character(RHS[[1]])

  ## store all parameters
  params <- character(length(parameter))
  for (i in seq_along(parameter)) {
    params[i] <- sprintf("PARAMETER(%s); \n", names(parameter)[i])
  }

  ## store all coefficents with linkfun
  if (!is.null(link)){
    coefs <- coefs_text <- character(length(link))
    for(i in seq_along(link)){
      coefs[i] <- trans_parnames(names(parameter)[i])
      ## filters out identity link
      if(!(coefs[i] %in% names(parameter))){
      coefs_text[i] <- sprintf("Type %s = %s; \n", coefs[i],
                          sprintf(all_links[[link[[i]]]], names(parameter)[i]))
      }
    }
  }

  ## variables without identity link
  trans_coefs <- setdiff(coefs, names(parameter))

  ## cpp script
  header <- "#include <TMB.hpp> \ntemplate<class Type> \nType objective_function<Type>::operator() () { \n"

  nll <- sprintf("Type nll = 0.0;\nnll = -sum(%s(%s, %s, %s, true));\nreturn nll;\n}",
                 as.character(RHS[[1]]), as.character(y), toString(RHS[2]),
                 paste(trans_coefs, collapse = ','))

  model <- paste0(header, paste(data_var, collapse=''),
                  paste(params, collapse=''),
                  paste(coefs_text, collapse=''), nll)

  ##return(model)
  write(model, file='template.cpp')
}

#' compiles template and create nll and gradient
#' @name TMB_mkfun
#' @param form A formula in expression form of "y ~ model"
#' @param parameter A list of initial values for the parameters
#' @param data A list of parameter in the formula with values in vectors
#' @param link link function and the model parameter
#' @importFrom TMB compile dynlib

TMB_mkfun <- function(formula,parameter,data=NULL,link=NULL){
  TMB_template(formula,parameter,data,link)
  TMB::compile("template.cpp")
  dyn.load(TMB::dynlib("template"))
  obj_fun <- MakeADFun(data=data, parameters=parameter, DLL="template")
  return(obj_fun)
}

