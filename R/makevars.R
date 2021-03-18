#' Make Makevars file to suppress irritating warnings
#' @name make_makevars
#' @examples
#' make_makevars()
#' @export
make_makevars <- function() {
  if (!dir.exists("~/.R")) {
    message("creating .R directory ...")
    dir.create("~/.R")
  }
  if (!file.exists("~/.R/Makevars")) {
    message("creating .R/Makevars ...")
    file.create("~/.R/Makevars")
  }
  message("adding CXXFLAGS to .R/Makevars ...")
  f <- file("~/.R/Makevars",open="a")
  writeLines(c("CXXFLAGS += -Wno-ignored-attributes"),
             con=f)
  close(f)
}
