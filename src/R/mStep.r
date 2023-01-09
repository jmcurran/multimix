
#' The M(aximisation) step
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
#'
#' @return an object of class \code{multimixParamList}---see 
#' \code{\link{initParamList}} for more information.
#' @export
mStep <- function(Z, D) {
  pistat <- colMeans(Z)
  W <- Z %*% diag(1/(D$n * pistat))  # Z scaled to have columns sum to 1 for use as weights.


  dstat <- list()
  for (i in seq_along(D$dlink)){ 
    dstat[[i]] <- crossprod(W, D$dvals[[i]])
  }

  ldstat <- list()
  for (i in seq_along(D$lcdisc)){ 
    ldstat[[i]] <- crossprod(W, D$ldvals[[i]])
  }

  ostat <- crossprod(W, D$ovals)
  ostat2 <- crossprod(W, D$ovals2)
  ovar <- ostat2 - ostat^2
  
  if(any(ovar <= .Machine$double.eps)){
    #browser()
    i = which(ovar <= .Machine$double.eps, arr.ind = TRUE)
    ovar[i] = .Machine$double.eps
    msg = paste0("Some variances have been replaced with machine precision to ",
                 "avoid divide by zero issues. This may have serious consequences.",
                 collapse = "\n")
    warning(msg)
  }

  cstat <- cstat2 <- cpstat <- list()
  for (i in seq_along(D$cdep)) {
    cstat[[i]] <- crossprod(W, D$cvals[[i]])
    cstat2[[i]] <- crossprod(W, D$cvals2[[i]])
    cpstat[[i]] <- crossprod(W, D$cprods[[i]])
  }

  lcstat <- lcstat2 <- lcpstat <- list()
  for (i in seq_along(D$lcdep)) {
    for (j in seq_len(D$ldlevs[i])) {
      group <- D$ldvals[[i]][, j] == 1
      gtot <- colSums(W[group, ])  ### or maybe pmin(colsums(W[group,]), minpstar)
      WW <- W[group, ] %*% diag(1/gtot)
      lcstat[[i]][[j]] <- crossprod(WW, D$lcvals[[i]][group, ])
      lcstat2[[i]][[j]] <- crossprod(WW, D$lcvals2[[i]][group, ])
      lcpstat[[i]][[j]] <- crossprod(WW, D$lcprods[[i]][group, ])
    }
  }

  cvar <- cstat
  ccov <- cpstat
  
  MVMV <- list()

  for (i in seq_along(D$cdep)) {
    lcdi <- length(D$cdep[[i]])
    nxp <- lcdi * (lcdi - 1)/2
    cvar[[i]] <- cstat2[[i]] - cstat[[i]]^2
    MVMV[[i]] <- list()

    for (j in seq_len(D$numClusters)) {
      MVMV[[i]][[j]] <- diag(length(D$cdep[[i]]))
      
      for (k in seq_along(D$cdep[[i]])) {
        MVMV[[i]][[j]][k, k] <- cvar[[i]][j, k]
      }
    }

    for (ii in seq_len(nxp)) {
      ccov[[i]][, ii] <- (cpstat[[i]][, ii] - cstat[[i]][, left(ii)] * cstat[[i]][, right(ii)])
    }

    for (j in seq_len(D$numClusters)) {
      for (ii in seq_len(nxp)) {
        MVMV[[i]][[j]][left(ii), right(ii)] <- MVMV[[i]][[j]][right(ii), left(ii)] <- ccov[[i]][j,
          ii]
      }
    }
  }

  lcvar <- lcstat
  lccov <- lcpstat
  LMV <- list()

  for (i in seq_along(D$lcdep)) {
    lcdi <- length(D$lcdep[[i]]) - 1
    nxp <- lcdi * (lcdi - 1)/2
    LMV[[i]] <- list()

    for (j in seq_len(D$numClusters)) {
      Temp <- rep(0, lcdi)
      LMV[[i]][[j]] <- diag(length(D$lcdep[[i]]) - 1)

      for (lv in seq_len(D$ldlevs[i])) {
        lcvar[[i]][[lv]] <- lcstat2[[i]][[lv]] - lcstat[[i]][[lv]]^2

        for (k in seq_len(lcdi)) {
          Temp[k] <- Temp[k] + lcvar[[i]][[lv]][j, k] * ldstat[[i]][j, lv]
        }  #k
      }  #lv  

      diag(LMV[[i]][[j]]) <- Temp
      M <- diag(Temp, nrow = length(Temp))
      for (lv in seq_len(D$ldlevs[i])) {
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
  
  P <- list(cpstat = cpstat,
            cstat = cstat,
            cstat2 = cstat2,
            cvar = cvar,
            dstat = dstat,
            lcpstat = lcpstat,
            lcstat = lcstat,
            lcstat2 = lcstat2,
            ldstat = ldstat,
            LMV = LMV,
            MVMV = MVMV,
            ostat = ostat,
            ovar = ovar,
            pistat = pistat,
            W = W
  )
  
  class(P) = "multimixParamList"
  return(P)
}
