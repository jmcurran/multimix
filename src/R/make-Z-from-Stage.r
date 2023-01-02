# Form an initial value of the Z matrix using the Stage 3/Stage 4 split

# Stage = scan(file = 'Stage.txt') Stage = Stage - 2 Init_grp = as.factor(Stage) Z <- model.matrix(~0 +
# Init_grp) attr(Z, 'assign') <- NULL attr(Z, 'contrasts') <- NULL colnames(Z) <- NULL rm(list = c('Stage',
# 'Init_grp')) W <- Z %*% diag(1/colSums(Z)) # Z scaled to have columns sum to 1 for use as weights.
