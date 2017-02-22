# ---- select-tree ----
trees <- read.tree(file = paste0(input.folder, 
                                 "tree/Stage2_MayrAll_Hackett_set10_decisive.tre"))
tree <- trees[[1]]
save(tree, file = paste0(input.folder, "tree/tree.RData"))