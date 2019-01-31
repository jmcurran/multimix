data_organise <- function(dframe, Init_grp, cdep=NULL, lcdep=NULL,
               minpstar = 1e-9, nIt = 60) {
Init_grp <- as.factor(Init_grp)              
## Functions 
count.unique <- function(x) length(unique(x))
# Pairing functions for small upper triangle / u < v /
pair.index <- function(u, v) {
   0.5*v^2 - 1.5*v + u + 1
}
right <- function(N) {
   floor(1.5 + 0.5*sqrt(-7 + 8*N))
}   
left <- function(N) {
   N - pair.index(1, right(N)) + 1
}
## Setup Constants 
n <- dim(dframe)[1]
p <- dim(dframe)[2]
q <- length(table(Init_grp))
if(length(Init_grp) != n) stop("initial grouping factor of wrong length")
discvar <- sapply(dframe, is.factor)  # logical indicator for discrete columns
#
ld <- rep(FALSE,p)
lc <- rep(FALSE,p)
for (i in seq_along(lcdep)) {
  ld[lcdep[[i]][1]] <- TRUE       # logical indicator for discrete columns in location model
  lc[lcdep[[i]][-1]] <- TRUE      # logical indicator for continuous columns in location model
}
md <- discvar & !ld
mc <- !discvar & !lc
oc <- mc
oc[unlist(cdep)] <- FALSE  # oc is dummy for continuous variables without local associations
olink <- seq_len(p)[oc]    # vector of continuous columns without local associations
dlink <- seq_len(p)[md]    # vector of discrete columns outside location cells
clink <- seq_len(p)[mc]    # vector of continuous columns outside location cells  
ldlink <- seq_len(p)[ld]   # vector of discrete columns inside location cells
lclink <- seq_len(p)[lc]   # vector of continuous columns inside location cells
lcdisc <- rep(NA, length(lcdep))
op <- length(olink)
for (i in seq_along(lcdep)) lcdisc[i] <- lcdep[[i]][1]  # discrete location columns in lcdep order
# number of levels of discrete variables:
dlevs <- apply(dframe[ , dlink, drop=FALSE], 2, count.unique) # needs all rows
ldlevs <- apply(dframe[ , lcdisc, drop=FALSE], 2, count.unique)
dvals <- list()
ldvals <- list()
for (i in seq_along(dlink) ) dvals[[i]] <- model.matrix(~ 0 + dframe[ ,dlink[i]])
for (i in seq_along(lcdisc) ) ldvals[[i]] <- model.matrix(~ 0 + dframe[ ,lcdisc[i]])
ovals <-  as.matrix(dframe[ ,olink])
ovals2 <-  ovals^2
cvals <- cvals2 <- list()
for (i in seq_along(cdep) ) {
  cvals[[i]] <- as.matrix(dframe[ ,cdep[[i]], drop=FALSE])
  cvals2[[i]] <- cvals[[i]]^2
}
lcvals <- lcvals2 <- list()
for (i in seq_along(lcdep) ) {
  lcvals[[i]] <- as.matrix(dframe[ ,lcdep[[i]][-1], drop=FALSE])
  lcvals2[[i]] <- lcvals[[i]]^2
}
cprods <- list()
for (i in seq_along(cdep) ) {
  lcdi <- length(cdep[[i]])
  nxp <- lcdi*(lcdi - 1)/2
  cprods[[i]] <- matrix(NA, nrow=n, ncol=nxp)
  for (ii in seq_len(nxp)) {
    cprods[[i]][ ,ii] <- as.matrix(dframe[ , cdep[[i]][left(ii)]]*dframe[ ,cdep[[i]][right(ii)], drop=FALSE])
  }}
lcprods <- list()
for (i in seq_along(lcdep) ) {
  lcdepi <- lcdep[[i]][-1]
  lcdi <- length(lcdepi)
  nxp <- lcdi*(lcdi - 1)/2
  lcprods[[i]] <- matrix(NA, nrow=n, ncol=nxp)
  for (ii in seq(length=nxp)) {
    lcprods[[i]][ ,ii]  <- 
      as.matrix(dframe[ , lcdepi[left(ii)]]*dframe[ ,lcdepi[right(ii)], drop=FALSE])
  }
}

ldxc <- list()
for (i in seq_along(lcdep) ) {
  ldxc[[i]] <- list()
  nlcd  <- max(0, length(lcdep[[i]]) - 1)
  for (j in seq_len(ldlevs[i])) {
    ldxc[[i]][[j]] <- matrix(NA, nrow=n, ncol=nlcd)
    for (k in seq_len(nlcd)) {
      ldxc[[i]][[j]][ , k] <- ldvals[[i]][ ,j]*lcvals[[i]][ , k]
    }}}
 D <- list(cdep = cdep,
           clink = clink,
           coltype = coltype,
           cprods = cprods,
           cvals = cvals,
           cvals2 = cvals2,
           dframe = dframe,
           discvar = discvar,
           dlevs = dlevs,
           dlink = dlink,
           dvals = dvals,
           Init_grp = Init_grp,
           lc = lc,
           lcdep = lcdep,
           lcdisc = lcdisc,
           lclink = lclink,
           lcprods = lcprods,
           lcvals = lcvals,
           lcvals2 = lcvals2,
           ld = ld,
           ldlevs = ldlevs,
           ldlink = ldlink,
           ldvals = ldvals,
           ldxc = ldxc,
           mc = mc,
           md = md,
           minpstar = minpstar,
           n = n,
           nIt = nIt,
           oc = oc,
           olink = olink,
           op = op,
           ovals = ovals,
           ovals2 = ovals2,
           p = p,
           q = q)
return(D)
}
