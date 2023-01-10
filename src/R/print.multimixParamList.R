#' S3 printing method for for multimix parameter results 
#'
#' @param x an object of class \code{multimixParamResults} -- see 
#' \code{\link{initParamList}} for more information.
#' @param type the statistic you want displayed. If \code{means} then the cluster
#' means will be displayed for each univariate continuous variable, the cluster proportions
#' for each level of a categorical variable, and the mean vector for each cluster and each
#' multivariate normal variable.
#' @param digits a vector of length 3. The first value determines how many decimal places to 
#' round categorical proportions to. The second value determines how many significant digits to 
#' display means to, and the third how many siginificant digits to display variances to. By default
#' proportions are rounded to 4 decimal places, means 2 significant digits, and variances 3 significant 
#' digits
#' @param ... additional arguments passed to \code{print}.
#'
#' @export
print.multimixParamList = function(x, type = c("means", "vars"), digits = c(4, 2, 3), ...){
  
  numCatVars = length(x$dstat)
  numUVNVars = ncol(x$ostat)
  numMVNVars = length(x$cstat)
  
  if(numCatVars > 0){
    cat(paste0("Categorical variables:\n"))
    cat(paste0("  ", paste0(names(x$dstat), collapse = ", "), "\n\n"))
  }
  
  if(numUVNVars > 0){
    cat(paste0("Univariate normal variables:\n"))
    cat(paste0("  ", paste0(colnames(x$ostat), collapse = ", "), "\n\n"))
  }
  
  if(numMVNVars > 0){
    cat(paste0("Multivariate normal variables:\n"))
    cells = sapply(x$cstat, function(cell){
                     paste0('(', paste0(colnames(cell), collapse = ", "), ')\n')
                  })
    cat(paste0("  Cell ", 1:length(cells), ": ", cells))
    cat("\n")
  }
  
  type = match.arg(type)
  
  if(type == "means"){
    if(numCatVars > 0){
      cat("Category proportions by cluster and level for categorical variables\n")
      cat("===================================================================\n")
      for(v in names(x$dstat)){
        cat(paste0("\nVariable ", v, ":\n"))
        print(round(x$dstat[[v]], digits = digits[1]))
        cat("\n")
      }
      cat("\n")
    }
    
    if(numUVNVars > 0){
      cat("Cluster means for univariate continuous variables\n")
      cat("===================================================================\n")
      for(v in colnames(x$ostat)){
        cat(paste0("\nVariable ", v, ":\n"))
        print(signif(x$ostat[,v, drop = FALSE], digits = digits[2]))
        cat("\n")
      }
      cat("\n")
    }
    
  }else{ ## vars
    
  }

}