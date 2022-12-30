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
#'   containing three elements: the \eqn{Z}{Z} matrix, the \eqn{P}{P} matrix,
#'   and a results matrix, called \code{results}, with \eqn{n}{n} rows and
#'   \eqn{numClusters}{numClusters} columns.
#' @export
#' @importFrom stats model.matrix runif
mmain <- function(D, Z, P, eps = 1e-9){
    nIter <- D$nIt
    numClusters <- ncol(Z)
    results <- matrix(0, nrow = nIter + 2, ncol = numClusters + 2)
    colnames(results)[c(1, numClusters + 2)] = c("Log-likelihood", "Iteration")
    colnames(results)[2:(numClusters + 1)] = paste0("Pr(G", 1:numClusters, ")")
    cyc <- 0
    zll <- P.to.Z(P, D)
    Z <- zll[[1]]
    llik <- zll[[2]]
    repeat {
        cyc <- cyc + 1
        P <- Z.to.P(Z, D, P)
        zll <- P.to.Z(P, D)
        Z <- zll[[1]]
        deltall <- zll[[2]] - llik
        if (deltall <= eps)
            break  # Explore choosing various epsilons, even 0.
        llik <- llik + deltall
        if (cyc >= nIter)
            break
        results[cyc + 1, ] <- c(llik, P$pistat, cyc)
    }  #cyc
    zpr <- list(Z = Z, P = P, results = results[1:cyc,])
    class(zpr) = "multimixResults"
    return(zpr)
}
