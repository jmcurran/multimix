#' Initialise the parameter list.
#'
#' Although the starting parameter list \code{P} may be specified directly, 
#. it is often tedious to do so.
#' Note also that any matrices specified must be positive definite.
#' This function calculates an initial \code{P} from \code{D} and a starting
#' value for \code{Z}.
#'
#' @param D an object of class \code{multimixSettings}---see
#'   \code{\link{data_organise}} for details.
#' @param Z an \eqn{n \times q}{n by q} matrix, where \eqn{n}{n} is the number
#'   of rows of \code{dframe} and \eqn{q}{q} is the number of components in the
#'   mixture. During the fitting \eqn{Z_{ij}}{Zij} holds the currently estimated
#'   probability that observation \eqn{i}{i} belongs to component \eqn{j}{j}.
#'   Often \code{Z} is initialized to a matrix of indicator columns for a
#'   partition of the data. It is also common to initialize \code{Z} to be the
#'   final \code{Z} from the fitting of a simpler model.
#'
#' @return an object of class \code{multimixParamList}---see 
#' \code{\link{Z.to.P}} for more information.
#' @export
initParamList <- function(D, Z) {
  len.cdep = length(D$cdep)
  dstat = vector("list", length(D$dlink))
  ldstat = vector("list", length(D$ldlink))
  cstat = vector("list", length(D$cdep))
  cstat2 = vector("list", len.cdep)
  cvar = vector("list", len.cdep)
  cpstat = vector("list", len.cdep)
  ccov = vector("list", len.cdep)


  lcstat <- lcstat2 <- lcpstat <- lcvar <- lccov <- replicate(length(D$lcdep), list())

  MVMV <- list()
  for (i in seq_along(D$cdep)) {
    MVMV[[i]] <- list()
    for (j in seq_len(D$numClusters)) {
      MVMV[[i]][[j]] <- diag(length(D$cdep[[i]]))
    }
  }

  LMV <- list()
  for (i in seq_along(D$lcdep)) {
    D$LMV[[i]] <- list()

    for (j in seq_len(D$numClusters)) {
      D$LMV[[i]][[j]] <- diag(length(D$lcdep[[i]]) - 1)
    }
  }

  pistat <- colSums(Z)/D$n 

  W <- Z %*% diag(1/(D$n * pistat)) # Matrix of weights

  for (i in seq_along(D$dlink)) {
    dstat[[i]] <- crossprod(W, D$dvals[[i]])
  }

  for (i in seq_along(D$lcdisc)) {
    ldstat[[i]] <- crossprod(W, D$ldvals[[i]])
  }

  ostat <- crossprod(W, D$ovals)
  ostat2 <- crossprod(W, D$ovals2)
  ovar <- ostat2 - ostat^2

  for (i in seq_along(D$cdep)) {
    cstat[[i]] <- crossprod(W, D$cvals[[i]])
    cstat2[[i]] <- crossprod(W, D$cvals2[[i]])
    cpstat[[i]] <- crossprod(W, D$cprods[[i]])
  }

  for (i in seq_along(D$lcdep)) {
    for (j in seq_len(D$ldlevs[i])) {
      group <- D$ldvals[[i]][, j] == 1
      gtot <- colSums(W[group, ])
      ### or maybe pmin(colsums(W[group,]), minpstar)
      WW <- W[group, ] %*% diag(1/gtot)
      lcstat[[i]][[j]] <- crossprod(WW, D$lcvals[[i]][group, ])
      lcstat2[[i]][[j]] <- crossprod(WW, D$lcvals2[[i]][group, ])
      lcpstat[[i]][[j]] <- crossprod(WW, D$lcprods[[i]][group, ])
    }
  }

  cvar <- cstat
  ccov <- cpstat

  for (i in seq_along(D$cdep)) {
    lcdi <- length(D$cdep[[i]])
    nxp <- lcdi * (lcdi - 1)/2
    cvar[[i]] <- cstat2[[i]] - cstat[[i]]^2

    for (j in 1:D$numClusters) {
      for (k in seq_along(D$cdep[[i]])) {
        MVMV[[i]][[j]][k, k] = cvar[[i]][j, k]
      }
    }

    for (ii in seq_len(nxp)) {
      ccov[[i]][, ii] <- (cpstat[[i]][, ii] - cstat[[i]][, left(ii)] * cstat[[i]][, right(ii)])
    }

    for (j in 1:D$numClusters) {
      for (ii in seq_len(nxp)) {
        MVMV[[i]][[j]][left(ii), right(ii)] <- MVMV[[i]][[j]][right(ii), left(ii)] <- ccov[[i]][j, ii]
      }
    }
  }

  lcvar <- lcstat
  lccov <- lcpstat

  for (i in seq_along(D$lcdep)) {
    lcdi <- length(D$lcdep[[i]]) - 1
    nxp <- lcdi * (lcdi - 1)/2

    for (j in 1:D$numClusters) {
      Temp <- rep(0, lcdi)

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
          M[right(ii), left(ii)] <- M[right(ii), left(ii)] + lccov[[i]][[lv]][j, ii] * ldstat[[i]][j, lv]
          M[left(ii), right(ii)] <- M[right(ii), left(ii)]
        }  #ii
      }  #lv
      LMV[[i]][[j]] <- M
    }  #j 
  }  #i

  P <- list(cpstat = cpstat, cstat = cstat, cstat2 = cstat2, cvar = cvar, dstat = dstat, lcpstat = lcpstat, lcstat = lcstat,
    lcstat2 = lcstat, ldstat = ldstat, LMV = LMV, MVMV = MVMV, ostat = ostat, ovar = ovar, pistat = pistat,
    W = W)

  class(P) = "multimixParamList"
  return(P)
}
