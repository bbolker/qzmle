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
#' set.seed(101)
#' rfp <- transform(emdbook::ReedfrogPred, nsize=as.numeric(size), random=rnorm(48))
#' form <- surv ~ dbinom(size = density,prob = exp(log_a)/(1 + exp(log_a)*h*density))
#' start <- list(h=1, log_a=c(2,0))
#' parameters <- list(log_a~poly(random))
#' links <- list(a="log")
#' TMB_template(form, start=start,links=links, parameters=parameters, data=rfp)


TMB_template <- function(formula,start,
                         links=NULL,
                         parameters=NULL,
                         data=NULL) {
  ## add data vectors
  if (!is.null(data)) {
    ## parse data in formula
    data_vars <- setdiff(all.vars(formula),
                         c(names(links), names(start), names(parameters)))

    data_vec <- character(length(data_vars))
    ## check data
    for (i in seq_along(data_vars)) {
        if (is.character(data$data_vars[i])) {
            stop("Cannot process string data")
        }
        ## store data vectors
        data_vec[i] <- sprintf("DATA_VECTOR(%s); \n", data_vars[i])
    }
  }

  ## parse distribution
  y <- formula[[2]]
  RHS <- RHS(formula)
  ddistn <- as.character(RHS[[1]])
  arglist <- as.list(RHS(formula)[-1])
  ## unname(sapply(arglist, deparse))

  ## submodels
  if(!missing(parameters)) {
    ## setting up submodels
    submodel_vars <- vapply(parameters,LHS_to_char,FUN.VALUE=character(1))
    Xlist <- lapply(parameters,parameter_parse, data=data)
    names(Xlist) <- submodel_vars
    ## make sure start values of parameters in the same order as the Xlist
    pvec <- start[submodel_vars]
  } else {
    submodel_vars <- NULL # starting values for submodel parameters
    pvec <- NULL # variables with submodel
  }

  ## if missing start arguments, use the named argument as the first element,
  ## all other elements of the sub-model parameter vector are 0
  for(i in submodel_vars){
    n_missed <- ncol(Xlist[[i]]) - length(pvec[[i]])
    if (n_missed < 0) stop('Too many argments in start for parameter: ', sQuote(i))
    if (n_missed != 0) pvec[[i]] <- c(pvec[[i]], rep(0, n_missed))
    ## add sub model parameter names
    names(pvec[[i]]) <- colnames(Xlist[[i]])
  }

  ## add parameters with no submodels
  start <- c(pvec, start[!names(start) %in% names(pvec)])

  ## store all parameters and submodel
  params <- character(length(start))
  submodel_eq <- submodel_data_text <- X_pname <- character(length(submodel_vars))

  for (i in seq_along(start)) {
    if (length(start[[i]]) > 1) {
      pname_param <- paste0(names(start)[i], "_param")
      params[i] <- sprintf("PARAMETER_VECTOR(%s); \n", pname_param[i])

      ## submodel data and parametrization
      if (names(start)[i] %in% submodel_vars) {
        X_pname <- paste0("X_", names(start)[i])
        submodel_data_text[i] <- sprintf("DATA_MATRIX(%s); \n", X_pname)
        submodel_eq[i] <- sprintf("vector <Type> %s = %s * %s",
                                  names(start)[i], X_pname[i], pname_param[i])
      }
    } else{
      params[i] <- sprintf("PARAMETER(%s); \n", names(start)[i])
    }
  }

  ## store all coefficents with linkfun
  ## vars names changed to link_varnames
  if (!is.null(links)){
    ## check if user place the same parameter as link function
    if(all(trans_parnames(names(start))==names(start))) {
      stop("parameter name in `start` should be in linkscale")

      ## FIXME: automatically change name?? necessary??
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

    coefs <- coefs_text <- character(length(links))
    trans_coefs <- character(length(links[links!="identity"]))
    for(i in seq_along(links)) {
      ## e.g "log_a" -> "a"
      coefs[i] <- trans_parnames(names(start)[i])
      ## filters out identity links (e.g. start has "log_a", a" is not in start)
      # links <- links[links!="identity"]
      if(!(coefs[i] %in% names(start))){
          coefs_len <- length(unlist(start[plinkfun(coefs[i], links[coefs[i]])]))
          link_pname <- sprintf(all_links[[links[[i]]]], names(start)[i])
          ## check if parameter with link is vector
          if (coefs_len > 1){
            coefs_text[i] <- sprintf("vector <Type> %s = %s; \n", coefs[i], link_pname)
          } else{
            coefs_text[i] <- sprintf("Type %s = %s; \n", coefs[i], link_pname)
          }
          trans_coefs[i] <- link_pname
          names(trans_coefs)[i] <- coefs[i]
        }
    }
  }


  ## make sure to vectorize parameter when needed
  for (i in seq_along(arglist)){
    if (length(grep(submodel_vars, arglist[[i]])) > 0) {
      arglist[[i]] <- sprintf("vector<Type> (%s)", deparse(arglist[[i]]))
    }
  }

  ## substitute transformed parameter without link name
  for (i in seq_along(trans_coefs)){
    arglist <- gsub(trans_coefs[[i]], names(trans_coefs)[i], arglist, fixed = T)
  }


  ## cpp script
  header <- "#include <TMB.hpp> \ntemplate<class Type> \nType objective_function<Type>::operator() () { \n"

  nll <- sprintf("Type nll = 0.0;\nnll -= sum(%s(%s, %s, true));\nreturn nll;\n}",
                 as.character(RHS[[1]]), as.character(y), paste(arglist, collapse = ','))

  model <- paste0(header,
                  paste(data_vec, collapse=''), ## data vectors
                  paste(submodel_data_text, collapse = ''), ## model matrix
                  paste(params, collapse=''), ## parameters
                  paste(submodel_eq, collapse = ''), ## submodel parametrization
                  paste(coefs_text, collapse=''), ## linkfun
                  nll)

  ##return(model)
  write(model, file='template.cpp')

  ## return list of data
  data_list <- list()
  for (i in data_vars){
    data_list[[i]] <- data[[i]]
  }
  names(Xlist) <- X_pname
  data_list <- c(data_list, Xlist)

  return(data=data_list)
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

  data_list <- TMB_template(formula,start,data,links)
  TMB::compile("template.cpp")
  dyn.load(TMB::dynlib("template"))
  obj_fun <- MakeADFun(data=data_list, parameters=start, DLL="template")
  return(obj_fun)
}

