#'  Read Z from Fortran output. Make into R matrix
#'
#' @param gr.out 
#'
#' @return
#' @export
#'
#' @examples
make_Z_fortran <- function(gr.out = "groups.out") {
  # Read Z from Fortran output. Make into R matrix.
  grout <- read.table(gr.out, header = FALSE)
    nc <- dim(grout)[2]
    Z <- as.matrix(grout[, -c(1, nc)])
    colnames(Z) <- NULL
    Z
}
