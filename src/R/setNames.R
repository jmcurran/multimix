#' @author James Curran
setPNames = function(P, D){
  numClusters = D$numClusters
  
  if(length(D$dlink) > 0){
    names(P$dstat) = names(D$dframe[, D$dlink])
    P$dstat = lapply(P$dstat, function(mat){
      colnames(mat) = gsub("^[^]]+[]]{2}(.*$)", "\\1", colnames(mat))
      rownames(mat) = paste0("C", 1:numClusters, ":")
      return(mat)
    })
  }
  
  if(length(D$olink) > 0){
    rownames(P$ostat) = paste0("C", 1:numClusters, ":")
    rownames(P$ovar) = paste0("C", 1:numClusters, ":")
  }
  
  for(cell in seq_along(D$cdep)){
    rownames(P$cstat[[cell]]) = paste0("C", 1:numClusters, ":")
    rownames(P$cvar[[cell]]) = paste0("C", 1:numClusters, ":")
    names(P$MVMV[[cell]]) = paste0("C", 1:numClusters, ":")
    for(cluster in seq_along(P$MVMV[[cell]])){
      rownames(P$MVMV[[cell]][[cluster]]) = 
        colnames(P$MVMV[[cell]][[cluster]]) = colnames(P$cstat[[cell]])
    }
  }
  
  if(length(D$lcdep) > 0){
    names(P$LMV) = names(P$lcstat) = names(D$dframe)[D$lcdisc]
    for(cell in seq_along(D$lcdep)){
      cellVarNames = names(D$dframe)[D$lcdep[[cell]]]
      discVarName = cellVarNames[1]
      contVarNames = cellVarNames[-1]
      discVarLevs = levels(D$dframe[, discVarName])
      
      names(P$LMV[[cell]]) = discVarLevs
      names(P$lcstat[[cell]]) = discVarLevs 
      
      for(lev in seq_along(P$LMV[[cell]])){
        rownames(P$LMV[[cell]][[lev]]) =
          colnames(P$LMV[[cell]][[lev]]) = contVarNames
        rownames(P$lcstat[[cell]][[lev]]) = paste0("C", 1:numClusters, ":")
      }
      
      # for(comp in seq_along(P$lcstat[[cell]])){
      #   rownames(P$lcstat[[cell]][[comp]]) = paste0("C", 1:numClusters, ":")
      # names(P$LMV[[cell]]) = paste0("C", 1:numClusters, ":")
      # for(cluster in seq_along(P$MVMV[[cell]])){
      #   rownames(P$MVMV[[cell]][[cluster]]) = 
      #     colnames(P$MVMV[[cell]][[cluster]]) = colnames(P$cstat[[cell]])
      #}
    }
  }
  
  return(P)
}
#' @author James Curran
setNames = function(mmResObj){
  numClusters = mmResObj$D$numClusters
  colnames(mmResObj$Z) = paste0("Pr(x[i,] in C", 1:numClusters, ")")
  
  mmResObj$P = setPNames(mmResObj$P, mmResObj$D)

  return(mmResObj)
}