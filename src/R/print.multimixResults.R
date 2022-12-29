#' S3 method for the printing of multimix results
#'
#' @param x an object of class \code{multimixResults}---see \code{\link{mmain}}
#'   for a description.
#' @param n display the last few iterations of the group probabilities. If 
#'   \code{TRUE} then the last 5 iterations will be displayed by default.
#'   Alternatively, a positive integer can be supplied. If this exceeds the 
#'   number of actual iterations, the output will be truncated.
#' @param ... other parameters passed to \code{print}. Not currently used.
#'
#' @importFrom methods is
#' @importFrom utils tail
#' @export
print.multimixResults = function(x, n = FALSE, ...){
  if(!is(x, "multimixResults")){
    stop("x must be of class multimixResults")
  }
  
  if(!is.logical(n) & !is.numeric(n)){
    stop("n must be either: TRUE, FALSE or an integer greater than zero.")
  }
  
  if(is.logical(n)){
    if(!n){
      numLines = 0
    }else{
      numLines = 5
    }
  }else{
    if(n <= 0){
      stop("n must be a positive integer")
    }else{
      numLines = floor(n)
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
  
  printProbs = function(probs){
    fmtString = paste0("%4d  ", 
                       paste0(rep("   %6.4f", numClusters), collapse = " "),
                       "\n")
    inputsString = paste0("probs[, ", 1:(numClusters + 1), "]", collapse = ", ")
    cmd = paste0("sprintf(\"", fmtString, "\", ", inputsString, ")")
    eval(parse(text = cmd))
  }
  
  if(numLines > 0){
    numLines = min(floor(numLines), nIterations)
    probRows = tail(x$results[, c(numClusters + 2, 2:(numClusters + 1))], numLines)
    cat(paste0("Last ", numLines, " values of the group probabilities:\n\n"))
    if(numClusters >= 10){
      cat(paste0("Iter.  ", paste0("Pr(G = ", sprintf("%2d", 1:numClusters), ")", collapse = " "), "\n"))
      cat(paste0("-----  ", paste0(rep("----------", numClusters), collapse = " "), "\n "))
    }else{
      cat(paste0("Iter.  ", paste0("Pr(G = ", 1:numClusters, ")", collapse = " "), "\n"))
      cat(paste0("-----  ", paste0(rep("---------", numClusters), collapse = " "), "\n "))
    }
    
    cat(printProbs(probRows))
  }
}