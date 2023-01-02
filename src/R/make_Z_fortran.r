#'  Read Z from FORTRAN output. Make into R matrix
#'
#' The FORTRAN version of Multimix produces two output files:
#' GENERAL.OUT and GROUPS.OUT. The latter mainly contains the Z matrix.
#'
#' This function facilitates the obtaining of Multimix R output given
#' Multimix FORTRAN output.
#'
#' @param gr.out string containing  a file name.
#'
#' @return a matrix containing a \eqn{Z}{Z} matrix.
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' Z <- make_Z_fortran(system.file('extdata', 'GROUPS-BP-Multimixf90.OUT', 
#'                     package = 'multimix'))
make_Z_fortran <- function(gr.out = "groups.out") {
  # Read Z from FORTRAN output. Make into R matrix.
  grout <- read.table(gr.out, header = FALSE)
  nc <- ncol(grout)
  Z <- as.matrix(grout[, -c(1, nc)])
  colnames(Z) <- NULL

  return(Z)
}

make_Z_fortrannew <- function(gr.out = "groups.out") {
  # Read Z from FORTRAN output. Make into R matrix.
  grout <- read.table(gr.out, header = FALSE)
  Init_grp = as.factor(grout[,1])
  Z <- model.matrix(~0 + Init_grp)
  attr(Z, "assign") <- NULL
  attr(Z, "contrasts") <- NULL
  colnames(Z) <- NULL
  
  return(Z)  
}

