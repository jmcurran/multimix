data("cancer.df")
D = data_organise(cancer.df, numClusters = 2)
stage = scan(system.file("extdata", "Stage.txt", package = "multimix")) - 2
Z <- make_Z_discrete(stage)
P <- initParamList(D, Z) 
zpr <- mmain(D, Z, P)
zpr
