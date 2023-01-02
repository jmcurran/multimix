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
count.unique <- function(x) {
  length(unique(x))
}

# Pairing functions for small upper triangle / u < v /
#' Maps integer pairs (u,v) with 0<u<v bijectively to positive integers.
#'
#' Used to reduce array dimensions by replacing
#'        A(x,y,z) by A*(x,pair.index(y,z))
#'
#' @param u positive integer scalar
#' @param v positive integer scalar
#'
#' @return integer scalar
#' @export
#'
#' @examples
#' pair.index(11,17)
#' pair.index(2,12)
pair.index <- function(u, v) {
  0.5 * v^2 - 1.5 * v + u + 1
}

#' Map integer index N>0 back to right member of generating pair.
#'
#' @param N  positive integer scalar
#'
#' @return  positive integer scalar
#' @export
#'
#' @examples
#' right(131)
#' right(57)
right <- function(N) {
  floor(1.5 + 0.5 * sqrt(-7 + 8 * N))
}


#' Map integer index N>0 back to left member of generating pair.
#'
#' @param N  positive integer scalar
#'
#' @return  positive integer scalar
#' @export
#'
#' @examples
#' left(131)
#' left(57)
left <- function(N) {
  N - pair.index(1, right(N)) + 1
}

