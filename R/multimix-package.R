#' multimix Model-based clustering using the EM (Expectation Maximisation) 
#' algorithm
#'
#' The package provides three categories of important functions:
#' \itemize{
#'  \item{operational}{ -- these functions are the functions used to perform 
#'  model fitting.}
#'  \item{helper}{ -- these functions are called internally and are unlikely to be
#'  called directly. They may not be exported in future versions of multimix.}
#'  \item{S3 methods}{ -- this set of functions helps with the display (printing
#'   or plotting) of the either the inputs or the results.}
#' }
#' 
#' @section multimix operational functions:
#' \code{data_organise}
#' \code{make_Z_discrete}
#' \code{make_Z_fortran}
#' \code{make_Z_random}
#' \code{initParamList}
#' \code{mmain}
#' \code{eStep}
#' \code{mStep}
#' 
#' @section multimix helper functions:
#' \code{count.unique}
#' \code{left}
#' \code{pair.index}
#' \code{right}
#' 
#' @section multmix S3 methods:
#' \code{plot.multimixResults}
#' \code{print.multimixParamList}
#' \code{print.multimixResults}
#'
#' @docType package
#' @name multimix
## usethis namespace: start
#' @importFrom Rcpp evalCpp
## usethis namespace: end
## usethis namespace: start
#' @useDynLib multimix, .registration = TRUE
## usethis namespace: end
NULL

NULL
