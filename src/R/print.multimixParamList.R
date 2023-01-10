#' S3 printing method for for multimix parameter results 
#'
#' @param x an object of class \code{multimixParamResults} -- see 
#' \code{\link{initParamList}} for more information.
#' @param type the statistic you want displayed. If \code{means} then the cluster
#' means will be displayed for each univariate continuous variable, the cluster proportions
#' for each level of a categorical variable, and the mean vector for each cluster and each
#' multivariate normal variable.
#' @param ... additional arguments passed to \code{print}.
#'
#' @export
print.multimixParamList = function(x, type = c("means", "vars"), ...){
  if(length(x$dstat) > 0){
    cat(paste0("Categorical variables:\n"))
    cat(paste0("  ", paste0(names(x$dstat), collapse = ", "), "\n"))
  }
  
  if(ncol(x$ostat) > 0){
    cat(paste0("\nUnivariate normal variables:\n"))
    cat(paste0("  ", paste0(colnames(x$ostat), collapse = ", "), "\n"))
  }
  
  if(length(x$cstat) > 0){
    cat(paste0("\nMultivariate normal variables:\n"))
    cells = sapply(x$cstat, function(cell){
                     paste0('(', paste0(colnames(cell), collapse = ", "), ')\n')
                  })
    cat(paste0("  Cell ", 1:length(cells), ": ", cells))
  }
  
  type = match.arg(type)
  
  if(type == "means"){
    
  }else{ ## vars
    
  }

}