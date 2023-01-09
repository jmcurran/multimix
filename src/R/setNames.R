setNames = function(mmResObj){
  numGroups = ncol(mmResObj$Z)
  colnames(mmResObj$Z) = paste0("Pr(obs(i) in G", 1:numGroups, ")")
  
  if(length(D$dlink) > 0){
    names(mmResObj$P$dstat) = names(D$dframe[, D$dlink])
    mmResObj$P$dstat = lapply(mmResObj$P$dstat, function(mat){
      colnames(mat) = gsub("^[^]]+[]]{2}(.*$)", "\\1", colnames(mat))
      return(mat)
    })
  }
  
  
  
  return(mmResObj)
}