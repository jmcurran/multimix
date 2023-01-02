#' Prepare data for use with multimix
#'
#' @param dframe a data frame containing the data set you wish to model.
#' @param numClusters the clusters you wish to fit.
#' @param numIter the maximum number of steps to that the EM agorithm will run
#' before terminating.
#' @param cdep a list of multivariate normal cells.
#' @param lcdep a list of location cells.
#' @param minpstar Minimum denominator for application of Bayes Rule.   
#'
#' @return An object of class \code{multimixSettings} which is a \code{list} 
#' with the following elements:
#' \itemize{
#'    \item{\code{cdep}}{a list of multivariate normal cells.}
#'    \item{\code{clink}}{column numbers of univariate normal variables.}
#'    \item{\code{cprods}}{a list over MVN cells containing a matrix of 
#'                         pair-wise products of columns in the cell, columns 
#'                         ordered by \code{\link{pair.index}}.}
#'    \item{\code{cvals}}{a list over MVN cells containing a matrix of columns of variables in the cell}
#'    \item{\code{cvals2}}{a list over MVN cells containing a matrix of squared columns of variables in the cell}
#'    \item{\code{dframe}}{the \code{data.frame} of variables}
#'    \item{\code{discvar}}{logical: the variable is takes values of either \code{TRUE} or \code{FALSE}}
#'    \item{\code{dlevs}}{for discrete cells: number of levels}
#'    \item{\code{dlink}}{column numbers of univariate discrete variables}
#'    \item{\code{dvals}}{a list over discrete cells of level indicator matrices}
#'    \item{\code{lc}}{logical: is continuous variable belonging to OT cell \code{TRUE}/\code{FALSE}}
#'    \item{\code{lcdep}}{a list of OT cells}
#'    \item{\code{lcdisc}}{column numbers of discrete variables in OT cells}
#'    \item{\code{lclink}}{column numbers of continuous variables in OT cells}
#'    \item{\code{lcprods}}{a list over OT cells containing a matrix of pair-wise products of continuous columns in the cell, columns ordered by \code{pair.index}}
#'    \item{\code{lcvals}}{a list over OT cells containing a matrix of continuous columns of variables in the cell}
#'    \item{\code{lcvals2}}{a list over OT cells containing a matrix of squared continuous columns of variables in the cell}
#'    \item{\code{ld}}{logical: is discrete variable belonging to OT cell \code{TRUE}/\code{FALSE}}
#'    \item{\code{ldlevs}}{for discrete variables in OT cells: number of levels}
#'    \item{\code{ldlink}}{column numbers of OT discrete variables}
#'    \item{\code{ldvals}}{a list over OT cells of level indicator matrices}
#'    \item{\code{ldxc}}{a list over OT cells whose members are lists over levels of matrices of the cell continuous variables whose columns are multiplied by the level indicator column}
#'    \item{\code{mc}}{logical: is continuous variable not in OT cell \code{TRUE}/\code{FALSE}}
#'    \item{\code{md}}{logical: is discrete variable not in OT cell \code{TRUE}/\code{FALSE}}
#'    \item{\code{minpstar}}{minimum denominator for appliction of Bayes' Rule}
#'    \item{\code{n}}{number of observations}
#'    \item{\code{numIter}}{the maximum number of steps to that the EM agorithm will run before terminating}
#'    \item{\code{oc}}{logical: is continuous variable in univariate cell \code{TRUE}/\code{FALSE}}
#'    \item{\code{olink}}{column numbers of continuous univariate cells}
#'    \item{\code{op}}{\code{length(olink)}}
#'    \item{\code{ovals}}{\code{n} by \code{op} matrix of continuous univariate variables}
#'    \item{\code{ovals2}}{\code{n} by \code{op} matrix of squared continuous univariate variables}
#'    \item{\code{numClusters}}{the number of clusters in the model.}
#' }
#' @export
#'
#' @examples
#' data(cancer.df)
#' mmObj = data_organise(cancer.df, numClusters = 2)
data_organise <- function(dframe, numClusters, numIter = 1000, cdep = NULL, lcdep = NULL, minpstar = 1e-09) {
  discvar <- sapply(dframe, is.factor)  # logical indicator for discrete columns
  #
  n <- nrow(dframe)
  v <- ncol(dframe)

  ld <- rep(FALSE, v)
  lc <- rep(FALSE, v)
  for (i in seq_along(lcdep)) {
    ld[lcdep[[i]][1]] <- TRUE  # logical indicator for discrete columns in location model
    lc[lcdep[[i]][-1]] <- TRUE  # logical indicator for continuous columns in location model
  }
  md <- discvar & !ld
  mc <- !discvar & !lc
  oc <- mc
  oc[unlist(cdep)] <- FALSE  # oc is dummy for continuous variables without local associations
  olink <- seq_len(v)[oc]  # vector of continuous columns without local associations
  dlink <- seq_len(v)[md]  # vector of discrete columns outside location cells
  clink <- seq_len(v)[mc]  # vector of continuous columns outside location cells  
  ldlink <- seq_len(v)[ld]  # vector of discrete columns inside location cells
  lclink <- seq_len(v)[lc]  # vector of continuous columns inside location cells
  lcdisc <- rep(NA, length(lcdep))
  op <- length(olink)
  for (i in seq_along(lcdep)) lcdisc[i] <- lcdep[[i]][1]  # discrete location columns in lcdep order
  # number of levels of discrete variables:
  dlevs <- apply(dframe[, dlink, drop = FALSE], 2, count.unique)  # needs all rows
  ldlevs <- apply(dframe[, lcdisc, drop = FALSE], 2, count.unique)
  dvals <- list()
  ldvals <- list()
  for (i in seq_along(dlink)) dvals[[i]] <- model.matrix(~0 + dframe[, dlink[i]])
  for (i in seq_along(lcdisc)) ldvals[[i]] <- model.matrix(~0 + dframe[, lcdisc[i]])
  ovals <- as.matrix(dframe[, olink])
  ovals2 <- ovals^2
  cvals <- cvals2 <- list()
  for (i in seq_along(cdep)) {
    cvals[[i]] <- as.matrix(dframe[, cdep[[i]], drop = FALSE])
    cvals2[[i]] <- cvals[[i]]^2
  }
  lcvals <- lcvals2 <- list()
  for (i in seq_along(lcdep)) {
    lcvals[[i]] <- as.matrix(dframe[, lcdep[[i]][-1], drop = FALSE])
    lcvals2[[i]] <- lcvals[[i]]^2
  }
  cprods <- list()
  for (i in seq_along(cdep)) {
    lcdi <- length(cdep[[i]])
    nxv <- lcdi * (lcdi - 1)/2
    cprods[[i]] <- matrix(NA, nrow = n, ncol = nxv)
    for (ii in seq_len(nxv)) {
      cprods[[i]][, ii] <- as.matrix(dframe[, cdep[[i]][left(ii)]] * dframe[, cdep[[i]][right(ii)], drop = FALSE])
    }
  }
  lcprods <- list()
  for (i in seq_along(lcdep)) {
    lcdepi <- lcdep[[i]][-1]
    lcdi <- length(lcdepi)
    nxv <- lcdi * (lcdi - 1)/2
    lcprods[[i]] <- matrix(NA, nrow = n, ncol = nxv)
    for (ii in seq(length = nxv)) {
      lcprods[[i]][, ii] <- as.matrix(dframe[, lcdepi[left(ii)]] * dframe[, lcdepi[right(ii)], drop = FALSE])
    }
  }

  ldxc <- list()
  for (i in seq_along(lcdep)) {
    ldxc[[i]] <- list()
    nlcd <- max(0, length(lcdep[[i]]) - 1)
    for (j in seq_len(ldlevs[i])) {
      ldxc[[i]][[j]] <- matrix(NA, nrow = n, ncol = nlcd)
      for (k in seq_len(nlcd)) {
        ldxc[[i]][[j]][, k] <- ldvals[[i]][, j] * lcvals[[i]][, k]
      }
    }
  }
  D <- list(cdep = cdep, clink = clink, cprods = cprods, cvals = cvals, cvals2 = cvals2, dframe = dframe, discvar = discvar,
    dlevs = dlevs, dlink = dlink, dvals = dvals, lc = lc, lcdep = lcdep, lcdisc = lcdisc, lclink = lclink,
    lcprods = lcprods, lcvals = lcvals, lcvals2 = lcvals2, ld = ld, ldlevs = ldlevs, ldlink = ldlink, ldvals = ldvals,
    ldxc = ldxc, mc = mc, md = md, minpstar = minpstar, n = n, numIter = numIter, oc = oc, olink = olink, op = op,
    ovals = ovals, ovals2 = ovals2, numClusters = numClusters)

  class(D) = "multimixSettings"
  return(D)
}
