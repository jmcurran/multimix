#' The M(aximisation) step
#'
#' @param P an object of class \code{multimixParamList}--see \code{\link{Z.to.P}} for more
#' information.
#' @param D an object of class \code{multimixSettings}---see
#'   \code{\link{data_organise}} for more information.
#'
#' @return a \code{list} containing two elements: a \code{matrix} named 
#' \code{Z}---see \code{\link{Z.to.P}} for more information, and a scalar
#' \code{llik} containing the current value of the log-likelihood.
#' @importFrom stats dnorm
#' @importFrom mvtnorm dmvnorm
#' @export
P.to.Z <- function(P, D) {
    with(c(P, D), {
      #browser()
        ollq <- matrix(0, nrow = n, ncol = numClusters)
        for (j in 1:numClusters) {
            ldens <- dnorm(as.vector(ovals), rep(ostat[j, ], rep(n, op)), rep(sqrt(ovar)[j, ], rep(n, op)),
                log = TRUE)
            lDENS <- matrix(ldens, nrow = n)
            ollq[, j] <- rowSums(lDENS)
        }

        cll <- matrix(0, nrow = n, ncol = numClusters)
        for (j in 1:numClusters) {
            cno <- 0
            while (cno < length(cdep)) {
                cno <- cno + 1
                cll[, j] <- cll[, j] + dmvnorm(cvals[[cno]], mean = cstat[[cno]][j, ], sigma = MVMV[[cno]][[j]],
                  log = TRUE)
            }
        }

        dllq <- matrix(0, nrow = n, ncol = numClusters)
        for (j in 1:numClusters) {
            for (v in seq_along(dstat)) {
                for (k in seq_len(dlevs[v])) {
                  dvk <- dvals[[v]][, k]
                  dllq[, j] <- dllq[, j] + dvk * log(dstat[[v]][j, k] + !dvk)
                }
            }
        }

        ldllq <- matrix(0, nrow = n, ncol = numClusters)
        for (j in 1:numClusters) {
            for (v in seq_along(ldstat)) {
                for (k in seq_len(ldlevs[v])) {
                  ldvk <- ldvals[[v]][, k]
                  ldllq[, j] <- ldllq[, j] + ldvk * log(ldstat[[v]][j, k] + !ldvk)
                }
            }
        }

        lcll <- matrix(0, nrow = n, ncol = numClusters)
        lmean <- list()
        for (j in 1:numClusters) {
            cno <- 0
            w <- W[, j]
            while (cno < length(lcdep)) {
                cno <- cno + 1
                nlev <- length(ldxc[[cno]])
                est <- vector("list", nlev)
                ndw <- colSums(diag(w) %*% ldvals[[cno]])
                for (lev in seq_len(nlev)) {
                  wdxc_ <- diag(w) %*% ldxc[[cno]][[lev]]
                  est[[lev]] <- colSums(wdxc_)/ndw[lev]
                }
                lmean[[cno]] <- ldvals[[cno]] %*% do.call(rbind, est)
                lcll[, j] <- lcll[, j] + dmvnorm(lcvals[[cno]] - lmean[[cno]], mean = rep(0, dim(lcvals[[cno]])[2]),
                  sigma = LMV[[cno]][[j]], log = TRUE)
            }
        }

        llx <- ollq + cll + dllq + ldllq + lcll
        rmx <- apply(llx, 1, max)
        expld <- exp(llx - rmx)
        pf_ <- t(expld) * pistat
        pstar <- colSums(pf_)
        pf_[, pstar < minpstar] <- 1/numClusters
        pstar[pstar < minpstar] <- 1
        Z <- t(pf_)/pstar
        llik <- sum(log(pstar) + rmx)
        zll <- list(Z = Z, llik = llik)
        return(zll)
    })
}

P.to.Znew <- function(P, D) {
  ollq <- 
    cll <- 
    dllq <- 
    ldllq <-
    lcll<-  matrix(0, nrow = D$n, ncol = D$numClusters)
  
  for (j in seq_len(D$numClusters)) {
    
    ## compute ollq
    ldens <- dnorm(as.vector(D$ovals), 
                   mean = rep(P$ostat[j, ], rep(D$n, D$op)), 
                   sd = rep(sqrt(P$ovar)[j, ], rep(D$n, D$op)),
                   log = TRUE)
    lDENS <- matrix(ldens, nrow = D$n)
    ollq[, j] <- rowSums(lDENS)

    ## compute CLL
    for (cno in seq_along(D$cdep)) {
      cll[, j] <- cll[, j] + dmvnorm(D$cvals[[cno]], 
                                     mean = P$cstat[[cno]][j, ], 
                                     sigma = P$MVMV[[cno]][[j]],
                                     log = TRUE)
    }
    
    
    ## computer dllq
    for (v in seq_along(dstat)) {
      for (k in seq_len(dlevs[v])) {
        dvk <- dvals[[v]][, k]
        dllq[, j] <- dllq[, j] + dvk * log(dstat[[v]][j, k] + !dvk)
      }
    }
    
    
    ## compute ldllq
    for (v in seq_along(ldstat)) {
      for (k in seq_len(ldlevs[v])) {
        ldvk <- ldvals[[v]][, k]
        ldllq[, j] <- ldllq[, j] + ldvk * log(ldstat[[v]][j, k] + !ldvk)
      }
    }
  }
  
  lmean <- list()
  for (j in seq_len(D$numClusters)) {
    cno <- 0
    w <- W[, j]
    while (cno < length(lcdep)) {
      cno <- cno + 1
      nlev <- length(ldxc[[cno]])
      est <- vector("list", nlev)
      ndw <- colSums(diag(w) %*% ldvals[[cno]])
      for (lev in seq_len(nlev)) {
        wdxc_ <- diag(w) %*% ldxc[[cno]][[lev]]
        est[[lev]] <- colSums(wdxc_)/ndw[lev]
      }
      lmean[[cno]] <- ldvals[[cno]] %*% do.call(rbind, est)
      lcll[, j] <- lcll[, j] + dmvnorm(lcvals[[cno]] - lmean[[cno]], mean = rep(0, dim(lcvals[[cno]])[2]),
                                       sigma = LMV[[cno]][[j]], log = TRUE)
    }
  }
  
  llx <- ollq + cll + dllq + ldllq + lcll
  rmx <- apply(llx, 1, max)
  expld <- exp(llx - rmx)
  pf_ <- t(expld) * pistat
  pstar <- colSums(pf_)
  pf_[, pstar < minpstar] <- 1/numClusters
  pstar[pstar < minpstar] <- 1
  Z <- t(pf_)/pstar
  llik <- sum(log(pstar) + rmx)
  zll <- list(Z = Z, llik = llik)

  return(zll)
}



