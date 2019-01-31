#' New Main Multimix Function
#'
#' @param D 
#'
#' @return
#' @export
#'
#' @examples
mmain <- function(D) {
  nIt <- D$nIt
  results <- matrix(0,nrow=D$nIt+1, ncol=D$q + 2)
  P <- first.Z.to.P(D)
  zll <- P.to.Z(P, D)
  Z <- zll[[1]]
  llik <- zll[[2]]
  results[1,] <- c(llik, P$pistat, 0)
  
  for (cyc in seq_len(nIt)){
    P <- Z.to.P(Z, D, P)
    zll <- P.to.Z(P, D)
    Z <- zll[[1]]
    llik <- zll[[2]]
    results[cyc+1,] <- c(llik, P$pistat, cyc)
  } #cyc
  
  zpr <- list(Z = Z, P = P, results = results)
  
  return(zpr)
}