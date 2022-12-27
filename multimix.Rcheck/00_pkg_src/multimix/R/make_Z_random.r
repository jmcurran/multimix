#' Start from random groups of similar size.
#'
#' A large number (\eqn{n}{n}) of observations are assigned randomly into
#' (\eqn{xq}{xq}) clusters. It is recommended to repeat Multimix runs with a
#' number of different seeds to search for a log-likelihood maximum.
#'
#' Also consider making additional clusters from observations with low
#' probabilities of belonging to any cluster in a previous clustering.
#'
#' @param xq   integer # Number of random "clusters" to be generated.
#' @param seed integer # Suggest using date as seed for random numbers.
#'
#' @return a matrix.
#' @importFrom runif
#' @export
#'
#' @examples
#' make_Z_random(3,271222)
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
	
    return(Z)
}
