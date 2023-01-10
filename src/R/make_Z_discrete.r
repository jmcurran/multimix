#' Make initial Z matrix from initial assignment of observations to clusters
#'
#' Z is an \eqn{n}{n} by \eqn{numClusters}{numClusters} matrix of non-negative numbers whose rows
#' sum to 1. The \eqn{ij^{\mathrm{th}}}{ij^th} element \eqn{z_{ij}}{z_ij} is a
#' probability that observation \eqn{i}{i} belongs to cluster \eqn{j}{j}. Rather
#' than begin from an initial assignment Multimix allows for a weighted
#' assignment accross several clusters.
#'
#' This function yields a 0/1 valued matrix.
#'
#' @param d integer
#'
#' @return a \code{matrix} whose entries are non-negative, and whose entries sum to 1.
#' @importFrom stats model.matrix
#' @export
#' @author Murray Jorgensen

#' @examples
#' stage = scan(file = system.file('extdata', 'Stage.txt', package = 'multimix'))
#' stage = stage - 2
#' Z = make_Z_discrete(stage)
make_Z_discrete <- function(d) {
  # Start from discrete variable with contiguous values lo:hi
  d1 <- d - min(d) + 1
  Init_grp = as.factor(d1)
  Z <- model.matrix(~0 + Init_grp)
  attr(Z, "assign") <- NULL
  attr(Z, "contrasts") <- NULL
  colnames(Z) <- NULL

  return(Z)
}
