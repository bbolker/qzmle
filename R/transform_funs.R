## 1. check names, warn if necessary
## (i.e. warn if (direction=="invlink" && name doesn't start with linkname)
## OR (direction=="link" && name starts with linkname)
## 2. transform x (with $invlink or $linkfun)
## 3. modify names(x) appropriately (add or subtract <linkname>_)

mklinkfun <- function(linkname, direction=c("linkfun","linkinv")) {
    direction <- match.arg(direction)
    m <- make.link(linkname)
    f <-function(x) {
      scale <- detect_scale(x)

      if (direction=="linkfun"){
        if (scale=="link"){
          warning("applying link functions: ",
                  paste(dQuote(linkname), "to possible link-scaled parameter:",
                        dQuote(names(x)), sep = " "))}
        x_name <- paste0(linkname,"_", names(x))
      }

      if (direction=="linkinv"){
        if (scale=="response"){
          warning("applying invlink functions: ",
                  paste(dQuote(linkname), "to possible not link-scaled parameter:",
                        dQuote(names(x)), sep = " "))}
        x_name <- trans_parnames(names(x))
      }

      result <- m[[direction]](x)
      names(result) <- x_name
      return(result)
    }
    return(f)
}

invlink <- mklinkfun("log")
x <- c(log_a=1)
invlink(x)


## checks whether name starts with either a specific link
## (or if NULL) starts with "some link"_
## returns either "response" or "link"
detect_scale <- function(x, link=NULL) {
  links <- names(x)
  regex <- sprintf("(%s)_", paste(names(all_links),collapse="|"))
  ifelse(length(grep(regex, links))>0,
         return("link"),
         return("response"))

  ## output a vector of link functions
  # s <- character(length(links))
  # for (i in seq_along(links)){
  #   ifelse(length(grep(regex, links[i]))>0,
  #          s[i]<- "link",
  #          s[i]<- "response")
  # }
  # return(s)
}

detect_scale(c(log_a=1)) # "link"
detect_scale(c(some_x=1)) # "response"
##detect_scale(c(log_a=1, some_x=1)) # "link" "response"




