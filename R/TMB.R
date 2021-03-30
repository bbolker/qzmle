#' Make model template to be compiled
#' @name TMB_template
#' @param formula A formula in expression form of "y ~ model"
#' @param start A list of initial values for the parameters
#' @param links link function and the model parameter
#' @param data A list of parameter in the formula with values in vectors
#' @examples
#' set.seed(123)
#' x <- runif(20, 1, 10)
#' y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
#' form <- y ~ dnorm(b0 + b1 * x, log_sigma)
#' links <- list(sigma="log")
#' start <- list(b0=1,b1=2, log_sigma=sd(y))
#' TMB_template(form, start=start,links=c(b0="identity", b1="identity", sigma="log"), data=list(x=x,y=y))
#' ff <- TMB_mkfun(form, start=start,links=c(b0="identity", b1="identity", sigma="log"), data=list(x=x,y=y))

TMB_template <- function(formula,start,links=NULL, parameters=NULL, data=NULL) {
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
  RHS <- RHS(formula)
  ddistn <- as.character(RHS[[1]])

  ## store all parameters
  params <- character(length(start))
  for (i in seq_along(start)) {
    params[i] <- sprintf("PARAMETER(%s); \n", names(start)[i])
  }

  ## store all coefficents with linkfun
  ## vars names changed to link_varnames
  if (!is.null(links)){
    ## check if user place the same parameter as link function
    if(all(trans_parnames(names(start))==names(start))) {
      stop("parameter name in `start` should be in linkscale")

      ## FIXME: automatically change name
      # link_ind <- which(names(start) %in% names(links))
      #  for(i in link_ind){
      #    ## change pnames to link_pnames using plinkfun
      #  }
      ## change start names
      # names(start) <- `[<-`(names(start), link_ind, )
    }

    ## if no links, add identity links
    if(length(links) != length(start)){
      pnames <- start[!trans_parnames(names(start)) %in% names(links)]
      ilinks <- rep(list("identity"), length(pnames))
      names(ilinks) <- names(pnames)
      links <- c(ilinks, links)
      links <- links[trans_parnames(names(start))]
    }

    coefs <- coefs_text <- character(sum(links!="identity"))
    for(i in seq_along(links)) {
      coefs[i] <- trans_parnames(names(start)[i])
      ## filters out identity links
      if(!(coefs[i] %in% names(start))){
        coefs_text[i] <- sprintf("Type %s = %s; \n", coefs[i],
                                 sprintf(all_links[[links[[i]]]], names(start)[i]))
        }
    # ## filters out identity links
    # links <- links[links!="identity"]
    # coefs <- link_coefs <- coefs_text <- character(length(links))
    # for(i in seq_along(links)) {
    #   coefs[i] <- names(links)[i]
    #   link_coefs[i] <- plinkfun(coefs[i], links[[i]])
    #   coefs_text[i] <- sprintf("Type %s = %s; \n", coefs[i],
    #                            sprintf(all_links[[links[[i]]]], link_coefs[i]))

    }
  }


  ## variables without identity links (unchanged var names)
  trans_coefs <- setdiff(coefs, names(start))

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


#' Compiles template and create nll and gradient
#' @name TMB_mkfun
#' @description compiles template and create nll and gradient
#' @param formula A formula in expression form of "y ~ model"
#' @param start A list of initial values for the parameters
#' @param data A list of parameter in the formula with values in vectors
#' @param links links function and the model parameter
#' @importFrom TMB compile dynlib MakeADFun

TMB_mkfun <- function(formula,start,links=NULL, parameters=NULL, data){

  TMB_template(formula,start,data,links)
  TMB::compile("template.cpp")
  dyn.load(TMB::dynlib("template"))
  obj_fun <- MakeADFun(data=data, parameters=start, DLL="template")
  return(obj_fun)
}

