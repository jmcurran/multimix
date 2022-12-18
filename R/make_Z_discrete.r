#' Title
#'
#' @param d 
#'
#' @return
#' @export
#'
#' @examples
make_Z_discrete <- function(d) {
  # Start from discrete variable with contiguous values lo:hi
  d1 <- d - min(d) + 1
    Init_grp = as.factor(d1)
    Z <- model.matrix(~0 + Init_grp)
    attr(Z, "assign") <- NULL
    attr(Z, "contrasts") <- NULL
    colnames(Z) <- NULL
    Z
}
