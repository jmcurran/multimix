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
#' @return an object of class \code{multimixParamList} which is a \code{list} 
#' with the following elements:
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
    LMV[[i]] <- list()

    for (j in seq_len(D$numClusters)) {
      LMV[[i]][[j]] <- diag(length(D$lcdep[[i]]) - 1)
    }
  }

  pistat <- colSums(Z) / D$n 

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
