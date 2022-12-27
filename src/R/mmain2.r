#' Title
#'
#' @param D an object of class \code{multimixSettings} - see
#'   \code{\link{data_organise}} for full description.
#' @param Z a matrix
#' @param P a matrix
#' @param eps
#'
#' @return an object of class \code{multimix results} which is a a list
#'   containing three elements: the \eqn{Z}{Z} matrix, the \eqn{P}{P} matrix,
#'   and a results matrix, called \code{results}, with \eqn{nIter + 2}{nIter +
#'   2} rows and \eqn{qq}{qq} columns.
#' @export
#' @importFrom stats model.matrix runif
#' @examples
mmain <- function(D, Z, P, eps = 1e-9){
    nIter <- D$nIt
    qq <- ncol(Z)
    results <- matrix(0, nrow = nIter + 2, ncol = qq + 2)
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
    zpr <- list(Z = Z, P = P, results = results)
    class(zpr) = "multimixResults"
    return(zpr)
}
