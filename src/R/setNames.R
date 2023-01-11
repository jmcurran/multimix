#' @author James Curran
setNames = function(mmResObj){
  numGroups = ncol(mmResObj$Z)
  colnames(mmResObj$Z) = paste0("Pr(x[i,] in G", 1:numGroups, ")")
  
  if(length(mmResObj$D$dlink) > 0){
    names(mmResObj$P$dstat) = names(mmResObj$D$dframe[, mmResObj$D$dlink])
    mmResObj$P$dstat = lapply(mmResObj$P$dstat, function(mat){
      colnames(mat) = gsub("^[^]]+[]]{2}(.*$)", "\\1", colnames(mat))
      rownames(mat) = paste0("G", 1:numGroups, ":")
      return(mat)
    })
  }
  
  if(length(mmResObj$D$olink) > 0){
    rownames(mmResObj$P$ostat) = paste0("G", 1:numGroups, ":")
    rownames(mmResObj$P$var) = paste0("G", 1:numGroups, ":")
  }
  
  for(cell in seq_along(mmResObj$P$cdep)){
    rownames(mmResObj$P$cstat[[cell]]) = paste0("G", 1:numGroups, ":")
    rownames(mmResObj$P$cvar[[cell]]) = paste0("G", 1:numGroups, ":")
    names(mmResObj$P$MVMV[[cell]]) = paste0("G", 1:numGroups, ":")
    for(group in seq_along(mmResObj$P$MVMV[[cell]])){
      rownames(mmResObj$P$MVMV[[cell]][[group]]) = 
        colnames(mmResObj$P$MVMV[[cell]][[group]]) = colnames(mmResObj$P$cstat[[cell]])
    }
  }
  
  if(length(mmResObj$D$lcdep) > 0){
    names(mmResObj$P$LMV) = names(mmResObj$D$dframe)[mmResObj$D$lcdisc]
    for(cell in seq_along(mmResObj$D$lcdep)){
      cellVarNames = names(D$dframe)[mmResObj$D$lcdep[[cell]]]
      discVarName = cellVarNames[1]
      contVarNames = cellVarNames[-1]
      discVarLevs = levels(mmResObj$D$dframe[, discVarName])
      
      names(mmResObj$P$LMV[[cell]]) = discVarLevs
      
      for(lev in seq_along(mmResObj$P$LMV[[cell]])){
        rownames(mmResObj$P$LMV[[cell]][[lev]]) =
          colnames(mmResObj$P$LMV[[cell]][[lev]]) = contVarNames
        rownames(mmResObj$P$lcstat[[cell]][[lev]]) = paste0("G", 1:numGroups, ":")
      }
      
      # for(comp in seq_along(mmResObj$P$lcstat[[cell]])){
      #   rownames(mmResObj$P$lcstat[[cell]][[comp]]) = paste0("G", 1:numGroups, ":")
        # names(mmResObj$P$LMV[[cell]]) = paste0("G", 1:numGroups, ":")
        # for(group in seq_along(mmResObj$P$MVMV[[cell]])){
        #   rownames(mmResObj$P$MVMV[[cell]][[group]]) = 
        #     colnames(mmResObj$P$MVMV[[cell]][[group]]) = colnames(mmResObj$P$cstat[[cell]])
      #}
    }
  }

  return(mmResObj)
}