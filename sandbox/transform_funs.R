function(linkname, direction=c("link","invlink")) {
    direction <- match.arg(direction)
    m <- make.link(linkname)
    f <-function(x) {
        ## 1. check names, warn if necessary
        ## (i.e. warn if (direction=="invlink" && name doesn't start with linkname) OR (direction=="link" && name starts with linkname)
        ## 2. transform x (with $invlink or $linkfun)
        ## 3. modify names(x) appropriately (add or subtract <linkname>_)
    }
    return(f)
}

detect_scale <- function(x, link=NULL) {
    ## checks whether name starts with either a specific link
    ## (or if NULL) starts with "some link"_
    ## returns either "response" or "link"
}

detect_scale(c(log_a=1)) # "link"
detect_scale(c(some_x=1)) # "response"

invlink <- mklinkfun("log")
invlink(x)



        
