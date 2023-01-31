#' S3 printing method for for multimix parameter results 
#'
#' @param x an object of class \code{multimixParamResults} -- see 
#' \code{\link{initParamList}} for more information.
#' @param type the statistic you want displayed. If \code{means} then the cluster
#' means will be displayed for each univariate continuous variable, the cluster proportions
#' for each level of a categorical variable, and the mean vector for each cluster and each
#' multivariate normal variable.
#' @param byLevel if \code{TRUE} then location model summary stats will be printed by the 
#' level of the factor in the location model. Otheriwse (default), they will be printed cluster
#' by cluster.
#' @param digits a vector of length 4. The first value determines how many decimal places to 
#' round categorical proportions to. The second value determines how many significant digits to 
#' display means to, and the third how many siginificant digits to display variances to. By default
#' proportions are rounded to 4 decimal places, means 2 significant digits, and variances 3 significant 
#' digits. The fourth value is only used if \code{pedantic == TRUE}, and is set to 16 significant figures by
#' default.
#' @param pedantic if \code{TRUE} then the results are printed to high precision for checking purposes.
#' This means \code{digits[4]} which is 16 decimal places by default.
#' @param raw if \code{TRUE} then switches off all of the customised printing and uses the default
#' print methods for \code{list}s etc.
#' @param ... additional arguments passed to \code{print}.
#' @returns No return value, called for side effects.
#' @author James Curran
#' @export
print.multimixParamList = function(x, type = c("means", "vars"), byLevel = FALSE, digits = c(4, 2, 3, 16), 
                                   pedantic = FALSE, raw = FALSE, ...){
  
  if(raw){
    print.AsIs(x)
  }else{
    numCatVars <- length(x$dstat)
    numUVNVars <- ncol(x$ostat)
    numMVNCells <- length(x$cstat)
    numLocationCells <- length(x$lcstat)
    
    if(numCatVars > 0){
      cat(paste0("Categorical variables:\n"))
      cat(paste0("  ", paste0(names(x$dstat), collapse = ", "), "\n\n"))
    }
    
    if(numUVNVars > 0){
      cat(paste0("Univariate normal variables:\n"))
      cat(paste0("  ", paste0(colnames(x$ostat), collapse = ", "), "\n\n"))
    }
    
    if(numMVNCells > 0){
      cat(paste0("Multivariate normal variables:\n"))
      cells <- sapply(x$cstat, function(cell){
                       paste0('(', paste0(colnames(cell), collapse = ", "), ')\n')
                    })
      cat(paste0("  Cell ", 1:length(cells), ": ", cells))
      cat("\n")
    }
    
    type <- match.arg(type)
    
    if(type == "means"){
      
      cat("")
      cat("Cluster proportions (pistat)\n")
      cat("============================\n")
      fmtString = if(!pedantic){
        paste0("%6.", digits[1], "f")
      }else{
        paste0("%.", digits[4], "g")
      }
      cat(paste0(paste0("C", 1:length(x$pistat), ": "), 
                 sprintf(fmtString, x$pistat), collapse = "\n"))
      cat("\n\n")
                 
      
      if(numCatVars > 0){
        cat("Category proportions by cluster and level for categorical variables\n")
        cat("===================================================================\n")
        
        for(v in names(x$dstat)){
          cat(paste0("\nVariable ", v, ":\n"))
          if(!pedantic){
            print(round(x$dstat[[v]], digits = digits[1]))
          }else{
            print(x$dstat[[v]], digits = digits[4])
          }
  #        cat("\n")
        }
        cat("\n")
      }
      
      if(numUVNVars > 0){
        cat("Cluster means for univariate continuous variables\n")
        cat("=================================================\n")
        for(v in colnames(x$ostat)){
          cat(paste0("\nVariable ", v, ":\n"))
          if(!pedantic){
            print(signif(x$ostat[,v, drop = FALSE], digits = digits[2]))
          }else{
            print(x$ostat[,v, drop = FALSE], digits = digits[4])
          }
  #        cat("\n")
        }
        cat("\n")
      }
      
      if(numMVNCells > 0){
        cat("Cluster means for multivariate continuous variables\n")
        cat("=================================================+=\n")
        for(cell in seq_along(x$cstat)){
          cat(paste0("\nCell ", cell, ":\n"))
          if(!pedantic){
            print(signif(x$cstat[[cell]], digits = digits[2]))
          }else{
            print(x$cstat[[cell]], digits = digits[4])
          }
  #        cat("\n")
        }
        cat("\n")
      }
      
      if(numLocationCells > 0){
        cat("Cluster means of  continuous variables of each component of the location model\n")
        cat("==============================================================================\n")
        
        cellNames = names(x$lcstat)
        if(byLevel){
          for(cell in seq_along(x$lcstat)){
            #browser()
            cat(paste0("\n", cellNames[cell], ":"))
            Levs = names(x$lcstat[[cell]])
            for(lev in seq_along(x$lcstat[[cell]])){
              cat(paste0("\n  ", Levs[lev], ":\n"))
              if(!pedantic){
                print(signif(x$lcstat[[cell]][[lev]], digits = digits[2]))
              }else{
                print(x$lcstat[[cell]][[lev]], digits = digits[4])
              }
  #            cat("\n")
            }
            cat("\n")
          }
        }else{
          for(cell in seq_along(x$lcstat)){
            Levs = names(x$lcstat[[cell]])
            numClusters = nrow(x$lcstat[[cell]][[1]]) # I don't like this but it will do for the time being
            #browser()
            
            for(cluster in seq_len(numClusters)){
              mat = do.call("rbind", lapply(x$lcstat[[cell]], function(m){m[cluster, , drop = FALSE]}))
              rownames(mat) = Levs
              cat(paste0("C", cluster, ":\n"))
              if(!pedantic){
                print(signif(mat, digits = digits[2]))
              }else{
                print(mat, digits = digits[4])
              }
  #            cat("\n")
            }
            cat("\n")
            #browser()
          }
        }
      }
      
    }else{ ## vars
      # if(numCatVars > 0){
      #   cat("Category proportions by cluster and level for categorical variables\n")
      #   cat("===================================================================\n")
      #   for(v in names(x$dstat)){
      #     cat(paste0("\nVariable ", v, ":\n"))
      #     print(round(x$dstat[[v]], digits = digits[1]))
      #     cat("\n")
      #   }
      #   cat("\n")
      # }
      
      if(numUVNVars > 0){
        cat("Cluster variances for univariate continuous variables\n")
        cat("================++++=================================\n")
        for(v in colnames(x$ovar)){
          cat(paste0("\nVariable ", v, ":\n"))
          if(!pedantic){
            print(signif(x$ovar[,v, drop = FALSE], digits = digits[2]))
          }else{
            print(x$ovar[,v, drop = FALSE], digits = digits[4])
          }
          cat("\n")
        }
        cat("\n")
      }
      
      if(numMVNCells > 0){
        cat("Cluster variances for multivariate continuous variables\n")
        cat("=======================================================\n")
        for(cell in seq_along(x$MVMV)){
          cat(paste0("\nCell ", cell, ":\n"))
          for(cluster in seq_along(x$MVMV[[cell]])){
            cat(paste0("Cluster ", cluster, ":\n"))
            if(!pedantic){
              print(signif(x$MVMV[[cell]][[cluster]], digits = digits[2]))
            }else{
              print(x$MVMV[[cell]][[cluster]], digits = digits[4])
            }
            cat("\n")
          }
        }
        cat("\n")
      }
      
      if(numLocationCells > 0){
        cat("Cluster variances of  continuous variables of each component of the location model\n")
        cat("==================================================================================\n")
        
        ##if(byLevel){
          cellNames = names(x$lcstat)
          for(cell in seq_along(x$lcstat)){
            #browser()
            cat(paste0("\n", cellNames[cell], ":"))
            Clusters = names(x$LMV[[cell]])
            for(cluster in seq_along(x$LMV[[cell]])){
              cat(paste0("\n  C", cluster, ":\n"))
              if(!pedantic){
                print(signif(x$LMV[[cell]][[cluster]], digits = digits[2]))
              }else{
                print(x$LMV[[cell]][[cluster]], digits = digits[4])
              }
              cat("\n")
            }
            cat("\n")
          }
        ##}else{ ## by cluster
          
        ##}
      }
    }
  }
}