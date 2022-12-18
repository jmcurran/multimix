## Helper functions

#' Count the number of unique items ion a vector x
#'
#' @param x a vector
#'
#' @return the number of unique items in \code{x}.
#' @export
#'
#' @examples
#' x = c(1, 2, 3)
#' count.unique(x)
#' 
#' x = c(1, 1, 1, 2, 3)
#' count.unique(x)
count.unique <- function(x){ 
  length(unique(x))
}

# Pairing functions for small upper triangle / u < v /
#' Title
#'
#' @param u 
#' @param v 
#'
#' @return
#' @export
#'
#' @examples
pair.index <- function(u, v) {
    0.5 * v^2 - 1.5 * v + u + 1
}

#' Title
#'
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
right <- function(N) {
    floor(1.5 + 0.5 * sqrt(-7 + 8 * N))
}


#' Title
#'
#' @param N 
#'
#' @return
#' @export
#'
#' @examples
left <- function(N) {
    N - pair.index(1, right(N)) + 1
}

