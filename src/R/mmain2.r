#' Title
#'
#' @param D 
#' @param Z 
#' @param P 
#'
#' @return
#' @export
#'
#' @examples
mmain <- function(D, Z, P) {
    results <- matrix(0, nrow = nIt + 2, ncol = qq + 2)
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
        if (cyc >= nIt)
            break
        results[cyc + 1, ] <- c(llik, P$pistat, cyc)
    }  #cyc
    zpr <- list(Z = Z, P = P, results = results)
    return(zpr)
}
