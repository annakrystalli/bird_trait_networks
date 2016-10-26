impute_varerror <- function(test.var, data) {
    
    vars <- ms_vars[ms_vars %in% names(data) & ms_vars != test.var]
    test.data <- data[apply(data, 1, function(x){!all(is.na(x))}), vars]

    
    g <- list(data = test.data, 
              meta = metadata, 
              mgm_types = mgm_types, meta_types = meta_types, 
              log.vars = log.vars, 
              spp.list = data.frame(species = rownames(test.data)),
              phylo.match = phylo.match,
              tree = drop.tip(tree, setdiff(tree$tip.label, rownames(test.data))))
    
    source(paste(script.folder, "mgm_dataprep.R", sep = ""), local = T)
    
    #' impute-mgm
    g <- impute_missForest(g)

    error <- data.frame(nrmse = g$missF_OOBerror["NRMSE"], pfc = g$missF_OOBerror["PFC"], 
                      test.var = test.var, data.size = prod(dim(g$data[names(g$missF_data)])), 
                      NAs = sum(is.na(g$data[names(g$missF_data)])), 
                      propNAs = sum(is.na(g$data[names(g$missF_data)]))/prod(dim(g$missF_data)),
                      nspp = nrow(g$missF_data), impt_var.no = ncol(g$missF_data))
  
    return(error)
}



#' ### impute missing data
#' mgm requires commplete cases or species are removed. Use `missForest`. 
#' Categorical and binary variables need to be supplied as vectors.
#+ impute-data, cache = T
impute_missForest <- function(g, cat.in = "numeric", cat.out = "numeric"){
  
  if(cat.in == "numeric"){  
    for(var in names(g$data)[g$mgm_types[names(g$data)] == "c"]){
      print(var)
      g$data[,var] <- as.factor(g$data[,var])
    }
  }
  
  impt <- missForest(g$data)
  g$missF_data <- impt$ximp
  g$missF_OOBerror <- impt$OOBerror
  
# return to numeric? Default = "numeric", ready for `mgmfit`.
  if(cat.out == "numeric"){
    for(var in names(g$missF_data)[g$mgm_types[names(g$missF_data)] == "c"]){
      g$missF_data[,var] <- as.numeric(g$missF_data[,var])
      g$data[,var] <- as.numeric(g$data[,var])
    }
  }
  
  return(g)
}

impute_phylopars <- function(g){
  
phy.data <- g$data[,mgm_types[gsub("log_", "", names(g$data))] != "c"]
phy.data <- phy.data[!apply(phy.data, 1, FUN = function(x){all(is.na(x))}),]
phy.data <- data.frame(species = rownames(phy.data), phy.data)

#13.38
t0 <- Sys.time()
PPE <- phylopars(trait_data = phy.data, 
                 tree = drop.tip(g$tree, setdiff(g$tree$tip.label, phy.data$species)),
                 pheno_error = F, pheno_correlated = F)

t1 <- Sys.time()
t1-t0
BIC(PPE)
}

