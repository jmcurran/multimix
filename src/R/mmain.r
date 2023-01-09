#' Title
#'
#' @param D an object of class \code{multimixSettings} - see
#'   \code{\link{data_organise}} for full description.
#' @param Z a matrix
#' @param P a matrix
#' @param eps Minimum increase in loglikelihood per EM step. If this is not
#'   exceeded the the algorithm will terminate.
#'
#' @return an object of class \code{multimix results} which is a a list
#'   containing four elements: the \code{multmixSettings} object \code{D},
#'   the \eqn{Z}{Z} matrix, the \eqn{P}{P} matrix,
#'   and a results matrix, called \code{results}, with \eqn{n}{n} rows and
#'   \eqn{numClusters}{numClusters} columns.
#' @export
#' @importFrom stats model.matrix runif
#' @examples
#' data(cancer.df)
#' D <- data_organise(cancer.df, numClusters = 2)
#' stage <- scan(system.file('extdata', 'Stage.txt', package = 'multimix')) - 2
#' Z <- make_Z_discrete(stage)
#' P <- initParamList(D,Z) 
#' zpr <- mmain(D,Z,P)
#' zpr
mmain <- function(D, Z, P, eps = 1e-09) {
  numIter <- D$numIter
  numClusters <- ncol(Z)
  results <- matrix(0, nrow = numIter + 2, ncol = numClusters + 2)
  colnames(results)[c(1, numClusters + 2)] = c("Log-likelihood", "Iteration")
  colnames(results)[2:(numClusters + 1)] = paste0("Pr(G", 1:numClusters, ")")
  cyc <- 0
  zll <- eStep(P, D)
  Z <- zll$Z
  llik <- zll$llik
  repeat {
    cyc <- cyc + 1

    P <- mStep(Z, D)
    zll <- eStep(P, D)

    #if(cyc == 31) browser()
    Z <- zll$Z
    #cat(paste0("Iteration: ", cyc, " LL: ", zll$llik, "\n"))
    deltall <- zll$llik - llik

    if (deltall <= eps)
      break  # Explore choosing various epsilons, even 0.

    llik <- llik + deltall

    if (cyc >= numIter)
      break

    results[cyc, ] <- c(llik, P$pistat, cyc)
  }  #cyc
  zpr <- setNames(list(D = D, Z = Z, P = P, results = results[1:(cyc - 1), ]))
  class(zpr) = "multimixResults"
  return(zpr)
}
