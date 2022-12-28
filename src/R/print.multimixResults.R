#' S3 method for the printing of multimix results
#'
#' @param x an object of class \code{multimixResults}---see \code{\link{mmain}}
#'   for a description.
#' @param dispLL display the last few iterations of the log-likelihood. If 
#'   \code{TRUE} then the last 5 iterations will be displayed by default.
#'   Alternatively, a positive integer can be supplied. If this exceeds the 
#'   number of actual iterations, the output will be truncated.
#' @param ... other parameters passed to \code{print}. Not currently used.
#'
#' @return 
#' @export
#'
#' @examples
print.multimixResults = function(x, dispLL = FALSE, ...){
  if(!is(x, "multimixResults")){
    stop("x must be of class multimixResults")
  }
  
  if(!is.logical(dispLL) & !is.numeric(dispLL)){
    stop("dispLL must be either: TRUE, FALSE or an integer greater than zero.")
  }
  
  if(is.logical(dispLL)){
    if(!dispLL){
      numLLlines = 0
    }else{
      numLLlines = 5
    }
  }else{
    if(dispLL <= 0){
      stop("dispLL must be a positive integer")
    }else{
      numLLlines = floor(dispLL)
    }
  }
  
  finalResult = x$results[nrow(x$results),]
  nIterations = finalResult["Iteration"]
  cat(paste0("The algorithm terminated after ", nIterations, 
             " iterations\n"))
  cat(paste0("The final value of the log-likelihood is ", 
            sprintf("%.2f", finalResult["Log-likelihood"]), "\n\n"))
  
  cat(paste0("The estimated group probabilities are:\n\n"))
  cat("i     Pr(G=i)\n")
  cat("----  -------\n")
  
  numClusters = length(finalResult) - 2
  for(i in 1:numClusters){
    cat(sprintf("%-4d  %6.4f\n", i, finalResult[i + 1]))
  }
  cat("\n")
  
  if(numLLlines > 0){
    numLines = min(floor(numLLlines), nIterations)
    llLines = tail(x$results[, c("Iteration", "Log-likelihood")], numLLlines)
    cat(paste0("Last ", numLines, " values of the log-likelihood:\n\n"))
    cat("Iter.  Log-lik.\n")
    cat("-----  --------\n")
    cat(paste0(sprintf("%5d %.2f", llLines[,1], llLines[,2]),
               collapse = "\n"))
  }
}