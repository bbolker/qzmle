#' Make model template to be compiled and return datalist
#' @name TMB_template
#' @param formula A formula in expression form of "y ~ model"
#' @param start A list of initial values for the parameters
#' @param links link function and the model parameter
#' @param parameters A list of linear submodels and random effects
#' @param data A list of parameter in the formula with values in vectors
#' @examples
#' set.seed(123)
#' x <- runif(20, 1, 10)
#' y <- rnorm(20, mean = 1.8 + 2.4 * x, sd = exp(0.3))
#' form <- y ~ dnorm(b0 + b1 * x, log_sigma)
#' links <- c(b0="identity", b1="identity", sigma="log")
#' start <- list(b0=1,b1=2, log_sigma=sd(y))
#' TMB_template(form, start=start,links=links, data=list(x=x,y=y))
#' ff <- TMB_mkfun(form, start=start,links=links, data=list(x=x,y=y))
#' ## submodel examples
#' set.seed(101)
#' rfp <- transform(emdbook::ReedfrogPred, nsize=as.numeric(size), random=rnorm(48))
#' form <- surv ~ dbinom(size = density,prob = exp(log_a)/(1 + exp(log_a)*h*density))
#' start <- list(h=1, log_a=c(2,0))
#' links <- list(a="log")
#' parameters <- list(log_a~poly(random))
#' TMB_template(form, start=start,links=links, parameters=parameters, data=rfp)
#' ## RE example with simulate data
#' rfpsim <- expand.grid(density=1:20,block=factor(1:20))
#' true_logit_a <- -1
#' true_log_h <- -1
#' true_logit_a_sd <- 0.3  ## log(0.3) = -1.20
#' set.seed(101)
#' logit_a_blk <- rnorm(20, mean=true_logit_a, sd=true_logit_a_sd)
#' a <- plogis(logit_a_blk[rfpsim$block])
#' prob <- a/(1 + a*exp(true_log_h)*rfpsim$density)
#' rfpsim$killed <- rbinom(nrow(rfpsim),size=rfpsim$density, prob=prob)
#' form <- killed ~ dbinom(size = density,prob = plogis(logit_a)/(1 + plogis(logit_a)*exp(log_h)*density))
#' parameters <- list(logit_a ~ 1 + (1|block))
#' links <- list(a="logit", h="log")
#' start <- list(logit_a=c(0), log_h=0)
#' TMB_template(form, start=start,links=links, parameters=parameters, data=rfpsim)
#' @importFrom lme4 nobars findbars mkReTrms
#' @importFrom Matrix t

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
    submodel_RHS <- sapply(parameters, onesided_formula)
    fixed_terms <- sapply(submodel_RHS, lme4::nobars)
    random_terms <- lapply(submodel_RHS, lme4::findbars)

    ## X matrix for fixed effects
    Xlist <- lapply(fixed_terms,parameter_parse, data=data)
    names(Xlist) <- submodel_vars
    ## make sure start values of parameters in the same order as the Xlist
    pvec <- start[submodel_vars]
  } else {
    submodel_vars <- NULL # starting values for submodel parameters
    pvec <- NULL # variables with submodel
  }


  ## Setting up random terms
  if (!all(sapply(random_terms, is.null))) {

    names(random_terms) <- submodel_vars
    random_terms <- unlist(random_terms[which(!is.null(random_terms))])

    ## only allow random intercept for now
    rand_int <- sapply(random_terms, "[[", 2)
    if (!all(rand_int==1)) stop("can only implement random intercept for now")

    ## set up random param
    Zlist <- re_rand_val <- vector(mode = "list", length = length(random_terms))
    re_param <- re_param_vec <- re_data <- re_eq <- nll_pen <- character(length(random_terms))
    re_sd <- numeric(length(random_terms))

    for (i in seq_along(random_terms)){
      Zlist[[i]] <- re_parameter_parse(list(random_terms[[i]]), data=data)
      Z_pname <- paste0("Z_", names(random_terms)[i])
      names(Zlist)[i] <- Z_pname
      re_data[i] <- sprintf("DATA_SPARSE_MATRIX(%s); \n", Z_pname)

      ## random sd in log scale
      re_sdname <- paste0(names(random_terms)[i], "_logsd")
      re_sd[[i]] <- 0 ## assume sd=1; log(sd) = 0
      names(re_sd)[i] <- re_sdname
      re_param[i] <- sprintf("PARAMETER(%s); \n", re_sdname)


      ## start value of 0 for random effects
      re_rand <- paste0(names(random_terms)[i], "_rand")
      re_rand_val[[i]] <- rep(0, length(levels(data[[random_terms[[1]][[3]]]])))
      names(re_rand_val)[i] <- re_rand
      re_param_vec[i] <- sprintf("PARAMETER_VECTOR(%s); \n", re_rand)

      ## add random effect on predictor
      re_eq[i] <- sprintf("%s += exp(%s) * (%s * %s);\n", names(random_terms)[i],
                          re_sdname, Z_pname, re_rand)

      ## penalization on nll
      nll_pen[i] <- sprintf("nll -= sum(dnorm(%s, Type(0), Type(1), true));\n", re_rand)
    }
  } else {re_data <- re_param <- re_eq <- re_param_vec <- re_rand <- re_sdname <- nll_pen <- Zlist <- NULL}

  ##FIXME: re_rand should be vector

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
  params <- start_pname <- character(length(start))
  submodel_eq <- submodel_data_text <- X_pname <- character(length(submodel_vars))

  for (i in seq_along(start)) {
    if ((length(start[[i]]) > 1) || (names(start)[i] %in% submodel_vars)) {
      pname_param <- paste0(names(start)[i], "_param")
      start_pname[i] <- pname_param
      params[i] <- sprintf("PARAMETER_VECTOR(%s); \n", pname_param)
    } else{
      start_pname[i] <- names(start)[i]
      params[i] <- sprintf("PARAMETER(%s); \n", names(start)[i])
    }

    ## submodel data and parametrization
    if (names(start)[i] %in% submodel_vars) {
      X_pname <- paste0("X_", names(start)[i])
      submodel_data_text[i] <- sprintf("DATA_MATRIX(%s); \n", X_pname)
      submodel_eq[i] <- sprintf("vector <Type> %s = %s * %s; \n",
                                names(start)[i], X_pname[i], pname_param[i])
    }
  }

  # if(is.null(links)){
  #   links <- character(length(start))
  # }

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
          if ((coefs_len > 1) || (names(start)[i] %in% submodel_vars)){
            coefs_text[i] <- sprintf("vector <Type> %s = %s; \n", coefs[i], link_pname)
          } else{
            coefs_text[i] <- sprintf("Type %s = %s; \n", coefs[i], link_pname)
          }
          trans_coefs[i] <- link_pname
          names(trans_coefs)[i] <- coefs[i]
          ## incase it is logit
          trans_coefs[i] <- gsub("invlogit", "plogis", trans_coefs[i])
        }
    }
  }


  ## make sure to vectorize parameter when needed
  for (i in seq_along(arglist)){
    if (length(grep(submodel_vars, arglist[[i]])) > 0) {
      arglist[[i]] <- sprintf("vector<Type>(%s)", deparse(arglist[[i]]))
    }
  }

  ## substitute transformed parameter without link name
  for (i in seq_along(trans_coefs)){
    arglist <- gsub(trans_coefs[[i]], names(trans_coefs)[i], arglist, fixed = T)
  }


  ## cpp script
  header <- "#include <TMB.hpp> \ntemplate<class Type> \nType objective_function<Type>::operator() () { \n"

  nll <- sprintf("Type nll = 0.0;\nnll -= sum(%s(%s, %s, true));\n",
                 as.character(RHS[[1]]), as.character(y), paste(arglist, collapse = ','))

  return_text <- "\nreturn nll;}\n"

  model <- paste0(header,
                  paste(data_vec, collapse=''), ## data vectors
                  paste(submodel_data_text, collapse = ''), ## model matrix
                  paste(re_data, collapse = ''), ## sparse matrix
                  paste(params, collapse=''), ## parameters
                  paste(re_param_vec, collapse = ''), ## random parameter
                  paste(re_param, collapse = ''), ## random sd
                  paste(submodel_eq, collapse = ''), ## submodel parametrization
                  paste(re_eq, collapse = ''), ## add r.e to predictor
                  paste(coefs_text, collapse=''), ## linkfun
                  nll, nll_pen, return_text)

  ##return(model)
  write(model, file='template.cpp')

  ## return list of data
  data_list <- list()
  for (i in data_vars){
    data_list[[i]] <- data[[i]]
  }
  names(start) <- start_pname
  start <- sapply(start, unname)
  int_n <- length(start)
  ##re_rand
  re_sd=list(rep(0, 20))
  start <- c(start, res_sdname=re_sd, re_rand=0)
  names(start) <- c(names(start)[1:int_n],re_rand, re_sdname)

  names(Xlist) <- X_pname
  data_list <- c(data_list, Xlist, Zlist, named_list(re_rand))

  return(list(data=data_list, start=start))
}


#' Compiles template and create nll and gradient
#' @name TMB_mkfun
#' @description compiles template and create nll and gradient
#' @param formula A formula in expression form of "y ~ model"
#' @param start A list of initial values for the parameters
#' @param data A list of parameter in the formula with values in vectors
#' @param links links function and the model parameter
#' @importFrom TMB compile dynlib MakeADFun

TMB_mkfun <- function(formula,start, links=NULL, parameters=NULL, data){
  data_list <- TMB_template(formula,start,links, parameters, data)
  TMB::compile("template.cpp")
  dyn.load(TMB::dynlib("template"))
  obj_fun <- MakeADFun(data=data_list$data, parameters= data_list$start, silent = T,
                       DLL="template",
                       random=data_list$re_rand
                       )
  obj_fun <- c(obj_fun, start = list(data_list$start))
  return(obj_fun)
}

