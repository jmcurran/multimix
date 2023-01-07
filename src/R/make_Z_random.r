#' Start from random groups of similar size.
#'
#' A large number (\eqn{n}{n}) of observations are assigned randomly into
#' (\eqn{xq}{xq}) clusters. It is recommended to repeat Multimix runs with a
#' number of different seeds to search for a log-likelihood maximum.
#'
#' Also consider making additional clusters from observations with low
#' probabilities of belonging to any cluster in a previous clustering.
#'
#' @param D   an object of class \code{multimixSettings} -- see 
#' \code{\link{data_organise}} for more information.
#' @param seed a positive integer to use as a random number seed.
#'
#' @return a matrix of dimension \eqn{n\times q}{n by q} where 
#' \eqn{n}{n} is the number of observations in \code{D$dframe}
#' and \eqn{q}{q} is the number of clusters in the model as specified
#' by \code{D$numClusters}.
#' 
#' @importFrom stats runif
#' @export
#'
#' @examples
#' data(cancer.df)
#' D = data_organise(cancer.df, numClusters = 2)
#' Z = make_Z_random(D)
#' table(Z)
make_Z_random <- function(D, seed = NULL) {
  if(!is(D, "multimixSettings")){
    stop("D must be an object of class multimixSettings.")
  }
  
  if(!is.null(seed)){
    if(seed <= 0){
      stop("Seed needs to be a non-zero positive integer")
    }else{
      set.seed(seed)
    }
  }
  
  n = nrow(D$dframe)
  numClusters = D$numClusters
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
