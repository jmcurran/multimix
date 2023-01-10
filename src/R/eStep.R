#' The E(xpectation) step
#'
#' @param P an object of class \code{multimixParamList}--see 
#' \code{\link{initParamList}} for more information.
#' @param D an object of class \code{multimixSettings}---see
#'   \code{\link{data_organise}} for more information.
#'
#' @return a \code{list} containing two elements: a \code{matrix} named 
#' \code{Z}---see \code{\link{mStep}} for more information, and a scalar
#' \code{llik} containing the current value of the log-likelihood.
#' @importFrom stats dnorm
#' @importFrom mvtnorm dmvnorm
#' @author Murray Jorgensen
#' @export
eStep <- function(P, D) {
  ollq <- cll <- dllq <- ldllq <- lcll <- matrix(0, nrow = D$n, ncol = D$numClusters)

  for (j in seq_len(D$numClusters)) {

    ## compute ollq
    ldens <- dnorm(as.vector(D$ovals), 
                   mean = rep(P$ostat[j, ], rep(D$n, D$op)), 
                   sd = rep(sqrt(P$ovar[j, ]), rep(D$n, D$op)), 
                   log = TRUE)
    lDENS <- matrix(ldens, nrow = D$n)
    ollq[, j] <- rowSums(lDENS)

    ## compute CLL
    for (cno in seq_along(D$cdep)) {
      cll[, j] <- cll[, j] + dmvnorm(D$cvals[[cno]], mean = P$cstat[[cno]][j, ], sigma = P$MVMV[[cno]][[j]],
        log = TRUE)
    }


    ## computer dllq
    for (v in seq_along(P$dstat)) {
      for (k in seq_len(D$dlevs[v])) {
        dvk <- D$dvals[[v]][, k]
        dllq[, j] <- dllq[, j] + dvk * log(P$dstat[[v]][j, k] + !dvk)
      }
    }


    ## compute ldllq
    for (v in seq_along(P$ldstat)) {
      for (k in seq_len(D$ldlevs[v])) {
        ldvk <- D$ldvals[[v]][, k]
        ldllq[, j] <- ldllq[, j] + ldvk * log(P$ldstat[[v]][j, k] + !ldvk)
      }
    }

    lmean <- list()
    w <- P$W[, j]

    for (cno in seq_along(D$lcdep)) {
      nlev <- length(D$ldxc[[cno]])
      est <- vector("list", nlev)
      ndw <- colSums(diag(w) %*% D$ldvals[[cno]])
      for (lev in seq_len(nlev)) {
        wdxc_ <- diag(w) %*% D$ldxc[[cno]][[lev]]
        est[[lev]] <- colSums(wdxc_)/ndw[lev]
      }
      lmean[[cno]] <- D$ldvals[[cno]] %*% do.call(rbind, est)
      lcll[, j] <- lcll[, j] + dmvnorm(D$lcvals[[cno]] - lmean[[cno]], mean = rep(0, ncol(D$lcvals[[cno]])),
        sigma = P$LMV[[cno]][[j]], log = TRUE)
    }
  }

  llx <- ollq + cll + dllq + ldllq + lcll
  rmx <- apply(llx, 1, max)
  expld <- exp(llx - rmx)
  pf_ <- t(expld) * P$pistat
  pstar <- colSums(pf_)
  pf_[, pstar < D$minpstar] <- 1/D$numClusters
  pstar[pstar < D$minpstar] <- 1
  Z <- t(pf_)/pstar
  llik <- sum(log(pstar) + rmx)
  zll <- list(Z = Z, llik = llik)

  return(zll)
}



