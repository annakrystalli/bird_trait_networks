#' ### impute missing data
#' mgm requires commplete cases or species are removed. Use `missForest`. 
#' Categorical and binary variables need to be supplied as vectors.
#+ impute-data, cache = T
for(var in names(g$data)[g$mgm_types[names(g$data)] == "c"]){
  g$data[,var] <- as.factor(g$data[,var])
}
impt <- missForest(g$data)
g$imp_data <- impt$ximp
g$OOBerror <- impt$OOBerror
#' return to numeric for `mgmfit`.
#+ return, cache = T
for(var in names(g$imp_data)[g$mgm_types[names(g$imp_data)] == "c"]){
  g$imp_data[,var] <- as.numeric(g$imp_data[,var])
  g$data[,var] <- as.numeric(g$data[,var])
}


PPE <- phylopars(trait_data = data.frame(species = rownames(g$data), g$data), tree = g$tree)
PPE_ou <- phylopars(trait_data = data.frame(species = rownames(g$data), g$data), 
                 tree = g$tree, model = "mvOU")

BIC(PPE_OU)