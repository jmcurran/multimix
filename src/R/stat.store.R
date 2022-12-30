#' Not currently exported as this is experimental
#' Murray - I have rewritten this to eliminate the use of with
#' and to make it clear which values are being copied from D, and which are
#' being redefined locally.
stat.store <- function(D) {
  # Set up store structures for sufficient statistics.
  # Values given not important.
  
  len.cdep = length(D$cdep)
  
  P = list(
    dstat <- vector("list", length(D$dlink)),
    ldstat <- vector("list", length(D$ldlink)),
    cstat <- vector("list", len.cdep),
    cstat2 <- vector("list", len.cdep),
    cvar <- vector("list", length(cdep)),
    cpstat <- vector("list", length(cdep)),
    ccov <- vector("list", length(cdep)),
    LMV <- list(),
    ostat = D$ostat,
    ostat2 = D$ostat2,
    ovar = D$ovar,
    MVMV <- list(),
    pistat <- D$pistat,
    W <- D$W
  )
  
  P$lcstat <- 
    P$lcstat2 <- 
    P$lcvar <- 
    P$lcpstat <- 
    P$lccov <-
    replicate(length(lcdep), list())
  
    
  for (i in seq_along(D$cdep) ) {
    MVMV[[i]] <- list()
    LMV[[i]] <- list()
    
    for (j in seq_len(D$qq)) {
      MVMV[[i]][[j]] <- diag(length(cdep[[i]]))
    }}
  }
  
  return(P)
}