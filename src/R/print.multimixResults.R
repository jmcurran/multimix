#' S3 method for the printing of multimix results
#'
#' @param x an object of class \code{multimixResults}---see \code{\link{mmain}}
#'   for a description.
#' @param ... other parameters passed to \code{print}. Not currently used.
#'
#' @return 
#' @export
#'
#' @examples
print.multimixResults = function(x, ...){
  if(!is(x, "multimixResults")){
    stop("x must be of class multimixResults")
  }
  
  finalResult = x$results[nrow(x$results),]
  cat(paste0("The algorithm terminated after ", finalResult["Iteration"], 
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
  cat("\n\n")
}