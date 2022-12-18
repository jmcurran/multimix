#' Start from random groups of similar size.
#'
#' @param xq 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
make_Z_random <- function(xq, seed = 310322) {
    set.seed(seed)
    x <- runif(n, 0, xq)
    x <- (n * x + 0.1)/(n + 0.2)
    y <- ceiling(x)
    Init_grp = as.factor(y)
    Z <- model.matrix(~0 + Init_grp)
    attr(Z, "assign") <- NULL
    attr(Z, "contrasts") <- NULL
    colnames(Z) <- NULL
    Z
}
