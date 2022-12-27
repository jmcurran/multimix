#' Title
#'
#' @param D
#' @param Z
#' @param P
#' @param nIter the maximum number of iteratins that multimix is allowed to run
#'   before it terminates
#' @param eps
#'
#' @return a list containing three elements: the \eqn{Z}{Z} matrix, the
#'   \eqn{P}{P} matrix, and a results matrix, called \code{results}, with
#'   \eqn{nIter + 2}{nIter + 2} rows and \eqn{qq}{qq} columns.
#' @export
#'
#' @examples
mmain <- function(D, Z, P, nIter = 100, eps = 1e-9) {
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
    return(zpr)
}
