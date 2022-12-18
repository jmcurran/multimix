#' Title
#'
#' @param P 
#' @param D 
#'
#' @return
#' @export
#'
#' @examples
P.to.Z <- function(P, D) {
    with(c(P, D), {
        ollq <- matrix(0, nrow = n, ncol = qq)
        for (j in 1:qq) {
            ldens <- dnorm(as.vector(ovals), rep(ostat[j, ], rep(n, op)), rep(sqrt(ovar)[j, ], rep(n, op)),
                log = TRUE)
            lDENS <- matrix(ldens, nrow = n)
            ollq[, j] <- rowSums(lDENS)
        }

        cll <- matrix(0, nrow = n, ncol = qq)
        for (j in 1:qq) {
            cno <- 0
            while (cno < length(cdep)) {
                cno <- cno + 1
                cll[, j] <- cll[, j] + dmvnorm(cvals[[cno]], mean = cstat[[cno]][j, ], sigma = MVMV[[cno]][[j]],
                  log = TRUE)
            }
        }

        dllq <- matrix(0, nrow = n, ncol = qq)
        for (j in 1:qq) {
            for (v in seq_along(dstat)) {
                for (k in seq_len(dlevs[v])) {
                  dvk <- dvals[[v]][, k]
                  dllq[, j] <- dllq[, j] + dvk * log(dstat[[v]][j, k] + !dvk)
                }
            }
        }

        ldllq <- matrix(0, nrow = n, ncol = qq)
        for (j in 1:qq) {
            for (v in seq_along(ldstat)) {
                for (k in seq_len(ldlevs[v])) {
                  ldvk <- ldvals[[v]][, k]
                  ldllq[, j] <- ldllq[, j] + ldvk * log(ldstat[[v]][j, k] + !ldvk)
                }
            }
        }

        lcll <- matrix(0, nrow = n, ncol = qq)
        lmean <- list()
        for (j in 1:qq) {
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
        pf_[, pstar < minpstar] <- 1/qq
        pstar[pstar < minpstar] <- 1
        Z <- t(pf_)/pstar
        llik <- sum(log(pstar) + rmx)
        zll <- list(Z = Z, llik = llik)
        return(zll)
    })
}





