#' Start from random groups of similar size.
#'
#' A large number (\eqn{n}{n}) of observations are assigned randomly into
#' (\eqn{xq}{xq}) clusters. It is recommended to repeat Multimix runs with a
#' number of different seeds to search for a log-likelihood maximum.
#'
#' Also consider making additional clusters from observations with low
#' probabilities of belonging to any cluster in a previous clustering.
#'
#' @param n integer - the number of observations.
#' @param numClusters   integer # Number of random 'clusters' to be generated.
#' @param seed integer # Suggest using date as seed for random numbers.
#'
#' @return a matrix.
#' @importFrom stats runif
#' @export
#'
#' @examples
#' Z = make_Z_random(300, 2, 271222)
#' table(Z)
make_Z_random <- function(n, numClusters, seed = 310322) {
  set.seed(seed)
  x <- runif(n, 0, numClusters)
  x <- (n * x + 0.1)/(n + 0.2)
  y <- ceiling(x)
  Init_grp = as.factor(y)
  Z <- model.matrix(~0 + Init_grp)
  attr(Z, "assign") <- NULL
  attr(Z, "contrasts") <- NULL
  colnames(Z) <- NULL

  return(Z)
}
