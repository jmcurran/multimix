#' First expectation step.
#'
#' @param D an object of class \code{multimixSettings}---see
#'   \code{\link{data_organise}} for details.
#' @param Z an \eqn{n \times q}{n by q} matrix, where \eqn{n}{n} is the number
#'   of rows of \code{dframe} and \eqn{q}{q} is the number of components in the
#'   mixture. During the fitting \eqn{Z_{ij}}{Zij} holds the currently estimated
#'   probability that observation \eqn{i}{i} belongs to component \eqn{j}{j}.
#'   Commonly \code{Z} is initialized to a matrix of indicator columns for a
#'   partition of the data.
#'
#' @return an object of class \code{multimixParamList}---see \code{\link{Z.to.P}} for more
#' information.
#' @export
first.Z.to.P <- function(D, Z) {
    with(D, {
        ## Z <- make_Z_random(qq) attr(Z,'assign') <- NULL attr(Z,'contrasts') <- NULL colnames(Z) <-
        ## NULL Set up store structures for sufficient statistics.  Values given not important.
        dstat <- vector("list", length(dlink))
        ldstat <- vector("list", length(ldlink))
        cstat <- vector("list", length(cdep))
        cstat2 <- vector("list", length(cdep))
        cvar <- vector("list", length(cdep))
        cpstat <- vector("list", length(cdep))
        ccov <- vector("list", length(cdep))
        ###
        lcstat <- lcstat2 <- lcvar <- lcpstat <- lccov <- replicate(length(lcdep), list())
        ###
        MVMV <- list()
        for (i in seq_along(cdep)) {
            MVMV[[i]] <- list()
            for (j in seq_len(qq)) {
                MVMV[[i]][[j]] <- diag(length(cdep[[i]]))
            }
        }
        LMV <- list()
        for (i in seq_along(lcdep)) {
            LMV[[i]] <- list()
            for (j in seq_len(qq)) {
                LMV[[i]][[j]] <- diag(length(lcdep[[i]]) - 1)
            }
        }
        results <- matrix(0, nrow = nIt + 1, ncol = qq + 2)
        ppi <- colSums(Z)/n
        W <- Z %*% diag(1/{
            n * ppi
        })

        for (i in seq_along(dlink)) dstat[[i]] <- crossprod(W, dvals[[i]])
        for (i in seq_along(lcdisc)) ldstat[[i]] <- crossprod(W, ldvals[[i]])
        ostat <- crossprod(W, ovals)
        ostat2 <- crossprod(W, ovals2)
        ovar <- ostat2 - ostat^2
        pistat <- ppi
        for (i in seq_along(cdep)) {
            cstat[[i]] <- crossprod(W, cvals[[i]])
            cstat2[[i]] <- crossprod(W, cvals2[[i]])
            cpstat[[i]] <- crossprod(W, cprods[[i]])
        }
        for (i in seq_along(lcdep)) {
            for (j in seq_len(ldlevs[i])) {
                group <- ldvals[[i]][, j] == 1
                gtot <- colSums(W[group, ])
                ### or maybe pmin(colsums(W[group,]), minpstar)
                WW <- W[group, ] %*% diag(1/gtot)
                lcstat[[i]][[j]] <- crossprod(WW, lcvals[[i]][group, ])
                lcstat2[[i]][[j]] <- crossprod(WW, lcvals2[[i]][group, ])
                lcpstat[[i]][[j]] <- crossprod(WW, lcprods[[i]][group, ])
            }
        }

        cvar <- cstat
        ccov <- cpstat
        for (i in seq_along(cdep)) {
            lcdi <- length(cdep[[i]])
            nxp <- lcdi * (lcdi - 1)/2
            cvar[[i]] <- cstat2[[i]] - {
                cstat[[i]]
            }^2
            for (j in 1:qq) {
                for (k in seq_along(cdep[[i]])) {
                  MVMV[[i]][[j]][k, k] = cvar[[i]][j, k]
                }
            }
            for (ii in seq_len(nxp)) {
                ccov[[i]][, ii] <- (cpstat[[i]][, ii] - cstat[[i]][, left(ii)] * cstat[[i]][, right(ii)])
            }
            for (j in 1:qq) {
                for (ii in seq_len(nxp)) {
                  MVMV[[i]][[j]][left(ii), right(ii)] <- MVMV[[i]][[j]][right(ii), left(ii)] <- ccov[[i]][j,
                    ii]
                }
            }
        }
        #
        lcvar <- lcstat
        lccov <- lcpstat
        for (i in seq_along(lcdep)) {
            lcdi <- length(lcdep[[i]]) - 1
            nxp <- lcdi * (lcdi - 1)/2
            for (j in 1:qq) {
                Temp <- rep(0, lcdi)
                for (lv in seq_len(ldlevs[i])) {
                  lcvar[[i]][[lv]] <- lcstat2[[i]][[lv]] - {
                    lcstat[[i]][[lv]]
                  }^2
                  for (k in seq_len(lcdi)) {
                    Temp[k] <- Temp[k] + lcvar[[i]][[lv]][j, k] * ldstat[[i]][j, lv]
                  }  #k
                }  #lv  
                diag(LMV[[i]][[j]]) <- Temp
                M <- diag(Temp, nrow = length(Temp))
                for (lv in seq_len(ldlevs[i])) {
                  for (ii in seq_len(nxp)) {
                    lccov[[i]][[lv]][, ii] <- (lcpstat[[i]][[lv]][, ii] - lcstat[[i]][[lv]][, left(ii)] * lcstat[[i]][[lv]][,
                      right(ii)])
                    M[right(ii), left(ii)] <- M[right(ii), left(ii)] + lccov[[i]][[lv]][j, ii] * ldstat[[i]][j,
                      lv]
                    M[left(ii), right(ii)] <- M[right(ii), left(ii)]
                  }  #ii
                }  #lv
                LMV[[i]][[j]] <- M
            }  #j 
        }  #i
        P <- list(dstat = dstat, ldstat = ldstat, ostat = ostat, ostat2 = ostat2, ovar = ovar, pistat = pistat,
            cstat = cstat, cstat2 = cstat2, cvar = cvar, cpstat = cpstat, lcstat = lcstat, lcstat2 = lcstat2,
            lcpstat = lcpstat, MVMV = MVMV, LMV = LMV, W = W)
        class(P) = "multimixParamList"
        return(P)
    })
}
