#' The E(xpectation) step
#'
#' Uses the current group membership to estimate the probabilities.
#'
#' @param Z an \eqn{n \times q}{n by q} matrix, where \eqn{n}{n} is the number
#'   of rows of \code{dframe} and \eqn{q}{q} is the number of components in the
#'   mixture. During the fitting \eqn{Z_{ij}}{Zij} holds the currently estimated
#'   probability that observation \eqn{i}{i} belongs to component \eqn{j}{j}.
#'   Commonly \code{Z} is initialized to a matrix of indicator columns for a
#'   partition of the data.
#' @param D an object of class \code{multimixSettings}---see
#'   \code{\link{data_organise}} for details.
#' @param P an object of class \code{multimixParamList} which is a \code{list} with the following elements:
#' \itemize{
#'    \item{\code{dstat}}{ -- \code{list} of matrices for each discrete variable
#'     not included in a location model. The matrix for each discrete variable 
#'     is made up of a column of length\eqn{q}{q} for each level (value) of the 
#'     variable giving the expected proportion of each level (column) for each 
#'     mixture component (row). Rows sum to 1.}
#'    \item{\code{ldstat}}{ -- \code{list} of matrices for each discrete 
#'    variable within a location model. The matrix for each discrete variable is
#'     made up of a column of length\eqn{q}{q} for each level (value) of the 
#'     variable giving the expected proportion of each level (column) for each 
#'     mixture component (row). Rows sum to 1. }
#'    \item{\code{ostat}}{ -- \code{matrix}  with a column for each continuous 
#'    variable outside any location mode whose\eqn{q}{q} rows give the current 
#'    estimated mean for each mixture component. }
#'    \item{\code{ostat2}}{ -- \code{matrix} with a column for each continuous 
#'    variable outside any location mode whose\eqn{q}{q}rows give the current 
#'    estimated mean square for each mixture component. }
#'    \item{\code{osvar}}{ -- \code{matrix} with a column for each continuous
#'     variable outside any location mode whose\eqn{q}{q} rows give the current
#'      estimated variance for each mixture component. }
#'    \item{\code{cstat}}{ -- \code{list} with a member for each nontrivial, 
#'    fully continuous, partition cell, that is not including discrete cells or 
#'    cells listed in lcdep, each member being a \code{matrix} with a column for each 
#'    continuous variable in that cell, whose\eqn{q}{q} rows give the current 
#'    estimated mean for each mixture component. }
#'    \item{\code{cstat2}}{ -- \code{list} with a member for each nontrivial, 
#'    fully continuous, partition cell, each member being a \code{matrix} with a column
#'     for each continuous variable in that cell, whose\eqn{q}{q} rows give the 
#'     current estimated mean square for each mixture component. }
#'    \item{\code{cvar}}{ -- \code{list} with a member for each nontrivial, 
#'    fully continuous, partition cell, each member being av\code{matrix} with a column
#'     for each continuous variable in that cell, whose\eqn{q}{q} rows give the 
#'     current estimated variance for each mixture component. }
#'    \item{\code{cpstat}}{ -- \code{list} with a member for each nontrivial, 
#'    fully continuous, partition cell, each member being the \code{matrix} with rows 
#'    for each of the\eqn{q}{q} mixture components and columns for each pair of 
#'    continuous variables in that cell, as ordered by \code{\link{pair.index}}. 
#'    The matrix 
#'    elements are the currently expected products of the variable pairs 
#'    arranged by component and pair. }
#'    \item{\code{ccov}}{ -- \code{list} with a member for each nontrivial, 
#'    fully continuous, partition cell, each member being the \code{matrix} with rows 
#'    for each of the \eqn{q}{q} mixture components and columns for each pair of 
#'    continuous variables in that cell, as ordered by pair.index. The matrix 
#'    elements are the currently expected covariances of the variable pairs 
#'    arranged by component and pair.}
#'    \item{\code{MVMV}}{ -- \code{list} with a member for each nontrivial, 
#'    fully continuous, partition cell, each member being a list with members 
#'    for each of the \eqn{q}{q} mixture components whose values are the 
#'    covariance matrix estimates for that cell and component. }
#'    \item{\code{lcstat}}{ -- \code{list} with a member for location partition 
#'    cell, each member being a \code{matrix} with a column for each continuous 
#'    variable in that cell, whose\eqn{q}{q} rows give the current estimated
#'     mean for each mixture component. }
#'    \item{\code{lcstat2}}{ -- \code{list} with a member for location partition
#'     cell, each member being a \code{matrix} with a column for each continuous 
#'     variable in that cell, whose \eqn{q}{q} rows give the current estimated 
#'     mean square for each mixture component. }
#'    \item{\code{lcpstat}}{ -- \code{list} with a member for each location 
#'    cell, each member being thev\code{matrix} with rows for each of the 
#'    \eqn{q}{q} mixture components and columns for each pair of continuous 
#'    variables in that cell, as ordered by \code{\link{pair.index}}. The matrix 
#'    elements are the currently expected products of the variable pairs 
#'    arranged by component and pair. }
#'    \item{\code{lccov}}{ -- \code{list} with a member for each location cell, 
#'    each member being the matrix with rows for each of the \eqn{q}{q} mixture 
#'    components and columns for each pair of continuous variables in that cell,
#'     as ordered by \code{\link{pair.index}}. The matrix elements are the 
#'     currently estimated covariances of the variable pairs arranged by 
#'     component and pair. }
#'    \item{\code{ldxcstat}}{ -- \code{list} with a member for each location 
#'    partition cell, each member being a list with a member for each level of 
#'    the cell's discrete variable that member being a matrix of mean values of 
#'    the continuous variables for each level-class combination. }
#' }
#'
#' @return an object of class \code{multimixParamList}---see above.
#' @export
Z.to.P <- function(Z, D, P) {
    with(c(D, P), {
        ppi <- colSums(Z)/n
        W <- Z %*% diag(1/{
            n * ppi
        })  # Z scaled to have columns sum to 1 for use as weights.

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
                gtot <- colSums(W[group, ])  ### or maybe pmin(colsums(W[group,]), minpstar)
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
            for (j in 1:numClusters) {
                for (k in seq_along(cdep[[i]])) {
                  MVMV[[i]][[j]][k, k] = cvar[[i]][j, k]
                }
            }
            for (ii in seq_len(nxp)) {
                ccov[[i]][, ii] <- (cpstat[[i]][, ii] - cstat[[i]][, left(ii)] * cstat[[i]][, right(ii)])
            }
            for (j in 1:numClusters) {
                for (ii in seq_len(nxp)) {
                  MVMV[[i]][[j]][left(ii), right(ii)] <- MVMV[[i]][[j]][right(ii), left(ii)] <- ccov[[i]][j, ii]
                }
            }
        }
        #
        lcvar <- lcstat
        lccov <- lcpstat
        for (i in seq_along(lcdep)) {
            lcdi <- length(lcdep[[i]]) - 1
            nxp <- lcdi * (lcdi - 1)/2
            for (j in 1:numClusters) {
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
            cstat = cstat, cstat2 = cstat2, cvar = cvar, cpstat = cpstat, lcstat = lcstat, lcstat2 = lcstat2, lcpstat = lcpstat,
            MVMV = MVMV, LMV = LMV, W = W)
        class(P) <- "multimixParamList"
        return(P)
    })
}


Z.to.Pnew <- function(Z, D, P) {
  P$pistat <- colMeans(Z)
  P$W <- Z %*% diag(1 / (D$n * P$pistat))  # Z scaled to have columns sum to 1 for use as weights.
  
  for (i in seq_along(D$dlink)) 
    P$dstat[[i]] <- crossprod(P$W, D$dvals[[i]])
  
  for (i in seq_along(D$lcdisc)) 
    P$ldstat[[i]] <- crossprod(P$W, D$ldvals[[i]])
  
  P$ostat <- crossprod(P$W, D$ovals)
  P$ostat2 <- crossprod(P$W, D$ovals2)
  P$ovar <- P$ostat2 - P$ostat^2

  for (i in seq_along(D$cdep)) {
    P$cstat[[i]] <- crossprod(P$W, D$cvals[[i]])
    P$cstat2[[i]] <- crossprod(P$W, D$cvals2[[i]])
    P$cpstat[[i]] <- crossprod(P$W, D$cprods[[i]])
  }
  
  for (i in seq_along(D$lcdep)) {
    for (j in seq_len(D$ldlevs[i])) {
      group <- D$ldvals[[i]][, j] == 1
      gtot <- colSums(W[group, ])  ### or maybe pmin(colsums(W[group,]), minpstar)
      WW <- W[group, ] %*% diag(1 / gtot)
      P$lcstat[[i]][[j]] <- crossprod(WW, D$lcvals[[i]][group, ])
      P$lcstat2[[i]][[j]] <- crossprod(WW, D$lcvals2[[i]][group, ])
      P$lcpstat[[i]][[j]] <- crossprod(WW, D$lcprods[[i]][group, ])
    }
  }
  
  P$cvar <- P$cstat
  ccov <- P$cpstat
  
  for (i in seq_along(D$cdep)) {
    lcdi <- length(D$cdep[[i]])
    nxp <- lcdi * (lcdi - 1) / 2
    P$cvar[[i]] <- P$cstat2[[i]] - P$cstat[[i]]^2
    
    for (j in seq_len(D$numClusters)) {
      for (k in seq_along(D$cdep[[i]])) {
        P$MVMV[[i]][[j]][k, k] = P$cvar[[i]][j, k]
      }
    }
    
    for (ii in seq_len(nxp)) {
      ccov[[i]][, ii] <- (P$cpstat[[i]][, ii] - 
                            P$cstat[[i]][, left(ii)] * 
                            P$cstat[[i]][, right(ii)])
    }
  
    for (j in seq_len(D$numClusters)) {
      for (ii in seq_len(nxp)) {
        P$MVMV[[i]][[j]][left(ii), right(ii)] <- 
          P$MVMV[[i]][[j]][right(ii), left(ii)] <- 
          ccov[[i]][j, ii]
      }
    }
  }
  
  P$lcvar <- P$lcstat
  lccov <- P$lcpstat
  
  for (i in seq_along(D$lcdep)) {
    lcdi <- length(D$lcdep[[i]]) - 1
    nxp <- lcdi * (lcdi - 1) / 2
    
    for (j in seq_len(D$numClusters)) {
      Temp <- rep(0, lcdi)
      
      for (lv in seq_len(D$ldlevs[i])) {
        lcvar[[i]][[lv]] <- lcstat2[[i]][[lv]] - lcstat[[i]][[lv]]^2
        
        for (k in seq_len(lcdi)) {
          Temp[k] <- Temp[k] + lcvar[[i]][[lv]][j, k] * ldstat[[i]][j, lv]
        }  #k
      }  #lv  
  
      diag(P$LMV[[i]][[j]]) <- Temp
      M <- diag(Temp, nrow = length(Temp))
      for (lv in seq_len(D$ldlevs[i])) {
        for (ii in seq_len(nxp)) {
          lccov[[i]][[lv]][, ii] <- (P$lcpstat[[i]][[lv]][, ii] - 
                                       P$lcstat[[i]][[lv]][, left(ii)] *
                                       P$lcstat[[i]][[lv]][, right(ii)])
          M[right(ii), left(ii)] <- M[right(ii), left(ii)] + 
                                      lccov[[i]][[lv]][j, ii] * 
                                      P$ldstat[[i]][j, lv]
          M[left(ii), right(ii)] <- M[right(ii), left(ii)]
        }  #ii
      }  #lv
      P$LMV[[i]][[j]] <- M
    }  #j 
  }  #i
  
  return(P)
}
