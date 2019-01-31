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

