#' This file contains functions we are testing. As a consequence, none of them 
#' are exported, and will need the ::: operator to use them
myfacto <- function(num) {
  ##dyn.load("facto.so")
  retvals <- .Fortran("facto",n = as.integer(num), answer = as.integer(1))
  return(retvals$answer)
}
