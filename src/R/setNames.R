setNames = function(mmResObj){
  numGroups = ncol(mmResObj$Z)
  colnames(mmResObj$Z) = paste0("Pr(x[i,] in G", 1:numGroups, ")")
  
  if(length(mmResObj$D$dlink) > 0){
    names(mmResObj$P$dstat) = names(D$dframe[, D$dlink])
    mmResObj$P$dstat = lapply(mmResObj$P$dstat, function(mat){
      colnames(mat) = gsub("^[^]]+[]]{2}(.*$)", "\\1", colnames(mat))
      rownames(mat) = paste0("G", 1:numGroups, ":")
      return(mat)
    })
  }
  
  if(ncol(mmResObj$P$ostat) > 0){
    rownames(mmResObj$P$ostat) = paste0("G", 1:numGroups, ":")
  }

  return(mmResObj)
}