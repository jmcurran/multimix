data("cancer.df")
D = data_organise(cancer.df, numClusters = 2)
stage = scan(system.file("extdata", "Stage.txt", package = "multimix")) - 2
system.time({
  for(i in 1:100){
    Z <- make_Z_discrete(stage)
    P <- initParamList(D,Z) 
    #P1 <- multimix:::initParamListNew(D,Z) 
    zpr <- mmain(D,Z,P)
  }
})
zpr1 <- mmain(D,Z,P1)
zpr
zpr1
