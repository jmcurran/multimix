#' Prepare data for use with multimix
#'
#' @param dframe a data frame containing the data set you wish to model.
#' @param numClusters the clusters you wish to fit.
#' @param niter the maximum number of steps to that the EM agorithm will run
#' before terminating.
#' @param cdep a list of multivariate normal cells.
#' @param lcdep a list of location cells.
#' @param minpstar Minimum denominator for appliction of Bayes Rule.   
#'
#' @return An object of class \code{multimixSettings} which is a \code{list} 
#' with the following elements:
#' \itemize{
#'    \item{\code{cdep}}{}
#'    \item{\code{clink}}{}
#'    \item{\code{cprods}}{}
#'    \item{\code{cvals}}{}
#'    \item{\code{cvals2}}{}
#'    \item{\code{dframe}}{}
#'    \item{\code{discvar}}{}
#'    \item{\code{dlevs}}{}
#'    \item{\code{dlink}}{}
#'    \item{\code{dvals}}{}
#'    \item{\code{lc}}{}
#'    \item{\code{lcdep}}{}
#'    \item{\code{lcdisc}}{}
#'    \item{\code{lclink}}{}
#'    \item{\code{lcprods}}{}
#'    \item{\code{lcvals}}{}
#'    \item{\code{lcvals2}}{}
#'    \item{\code{ld}}{}
#'    \item{\code{ldlevs}}{}
#'    \item{\code{ldlink}}{}
#'    \item{\code{ldvals}}{}
#'    \item{\code{ldxc}}{}
#'    \item{\code{mc}}{}
#'    \item{\code{md}}{}
#'    \item{\code{minpstar}}{}
#'    \item{\code{n}}{}
#'    \item{\code{nIt}}{}
#'    \item{\code{oc}}{}
#'    \item{\code{olink}}{}
#'    \item{\code{op}}{}
#'    \item{\code{ovals}}{}
#'    \item{\code{ovals2}}{}
#'    \item{\code{qq}}{the number of clusters in the model.}
#' }
#' @export
#'
#' @examples
#' data(cancer.df)
#' mmObj = data_organise(cancer.df)
data_organise <- function(dframe,
                          numClusters,
                          niter  = 1000,
                          cdep = NULL, 
                          lcdep = NULL, 
                          minpstar = 1e-09) {
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
    D <- list(cdep = cdep, clink = clink, cprods = cprods, cvals = cvals, cvals2 = cvals2, dframe = dframe,
        discvar = discvar, dlevs = dlevs, dlink = dlink, dvals = dvals, lc = lc, lcdep = lcdep, lcdisc = lcdisc,
        lclink = lclink, lcprods = lcprods, lcvals = lcvals, lcvals2 = lcvals2, ld = ld, ldlevs = ldlevs, ldlink = ldlink,
        ldvals = ldvals, ldxc = ldxc, mc = mc, md = md, minpstar = minpstar, n = n, nIt = niter, oc = oc, olink = olink,
        op = op, ovals = ovals, ovals2 = ovals2, qq = numClusters)
    
    class(D) = "multimixSettings"
    return(D)
}
