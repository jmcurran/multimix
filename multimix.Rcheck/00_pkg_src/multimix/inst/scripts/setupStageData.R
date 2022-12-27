stage = scan(file = system.file("extdata", "Stage.txt", package = "multimix"))
stage = stage - 2
save(stage, file = "data/stage.Rda")
