# ---- select-tree ----
trees <- read.tree(file = "tree/Stage2_MayrAll_Hackett_set10_decisive.tre")
tree <- trees[[1]]
save(tree, file = "tree/tree.RData")