#' calculate number of variable pairwise data-point availability 
#'
#' @param data a wide data_frame containing variables in columns and species in rows. Function ignores columns named "species" and "synonyms". 
#' @details 
#' 
#' @return a data.frame of the unique variable pairwise combinations and counts $(n) of data availability. The dataframe is sorted according to $n highest to lowest
#' @export
#'
#' @examples
#' 
# ---- calcTraitPairN-f ----
calcTraitPairN <- function(data){
  
  vars <- names(data)
  
  var.grid <- expand.grid(vars, vars, stringsAsFactors = F)
  var.grid <- var.grid[var.grid[,1] != var.grid[,2],]
  
  indx <- !duplicated(t(apply(var.grid, 1, sort))) # finds non - duplicates in sorted rows
  var.grid <- var.grid[indx, ]
  
  countN <- function(x, data){sum(complete.cases(data[,c(x[1], x[2])]))}
  
  var.grid <- data.frame(var.grid, n = apply(var.grid, 1, FUN = countN, data = data))
  
}

pglsPhyloCor <- function(x, data, tree, log.vars = NULL, datTypes, 
                         result = "row"){
  # ---- log-log.vars ----  
  if(datTypes == "nn") {
    for(var in x){
      if(validateLog(var = var, log.vars, data)){
        data[, var] <- log(data[, var])
        setNames(data[, var], paste("log", var, sep = "_"))
        x[x == var] <- paste("log", var, sep = "_")
      }}}
  
  # ---- set-up-vars ----
  var1 <- unlist(x[1])
  var2 <- unlist(x[2])
  data <- data[, c("species", "synonyms", var1, var2)]
  
  # ---- get-TD ----  
  TD <- getTD(data, tree)
  
  # ---- sub-complete-cases ----
  data <- data[complete.cases(data),]
  
  # ---- fit-cor ---- 
  if(datTypes == "nn") {
    row <- fitNumDat(data, tree, var1, var2, TD, result)
  }
  if(datTypes == "bb") {
    row <-  fitBinDat(data, tree, var1, var2, TD, result) 
  }
  if(datTypes == "cc") {
    meta <- list(getLevels(metadata, var1), getLevels(metadata, var2))
    # factorise data
    data[,c(var1, var2)] <- mapply(FUN = function(x, meta){factor(as.character(x),
                                                                  levels = meta$levels,
                                                                  labels = meta$labels)},
                                   x = data[,c(var1, var2)], meta = meta) %>%
      as.data.frame(stringsAsFactors = T)
    
    row <-  fitCatDat(data, tree, var1, var2, TD, result) }
  
  # ---- return ----
  return(row)
  
}

fitNumDat <- function(data, tree, var1, var2, TD, result = "row") {
  
  # ---- configure-tree ----
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$synonyms))
  tree$tip.label <- data$species[match(tree$tip.label, data$synonyms)] 
  
  # ---- get-phylosig ----
  physig.var1 <- phylosig(tree, setNames(data[,var1], data$species))
  physig.var2 <- phylosig(tree, setNames(data[,var2], data$species))
  
  
  # ---- get-cor ----
  cor <- cor(data[,var1], data[,var2])

  ## ---- get-model ----
  cd <- comparative.data(phy = tree, data = data, names.col = "synonyms", vcv=F)
  
  result.pgls <- try(pgls(as.formula(paste(var1, "~", var2, sep = "")),
                          data = cd, lambda="ML"))
  
  # ---- return-model ----
  if(result == "model"){
    return(results.pgls)
  }
  
  # ---- return-row ---- 
  if(class(result.pgls) == "try-error"){
    phylocor2 <- NA
    lambda <- NA
    p <- NA
    ci <- NA
    error <- gsub(pattern = ".*\n  ", "", geterrmessage())
    }else{
    t <- summary(result.pgls)$coefficients[var2,3]
    df <- as.vector(summary(result.pgls)$fstatistic["dendf"])
    phylocor2 <- sqrt((t*t)/((t*t)+df))*sign(summary(result.pgls)$coefficients[var2,1])
    lambda <- result.pgls$param["lambda"]
    p <- coef(summary(result.pgls))[2,4]
    ci <- result.pgls$param.CI$lambda$ci
    error <- NA
  }
  return(data.frame(var1 = var1, var2 = var2, n = length(data$species),
                    Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                    physig.var1, physig.var2,
                    cor = cor, phylocor = phylocor2, 
                    lambda = lambda, p = p, ci = ci, 
                    error = error))
}

fitBinDat <- function(data, tree, var1, var2, TD, result = "row") {

  # ---- configure-tree ----
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$synonyms))
  tree$tip.label <- data$species[match(tree$tip.label, data$synonyms)]  
  
  # ---- get-cor ----
  cor <- NA
  
  # ---- get-cor ----
  data <- data[match(tree$tip.label, data$species),]
  x <- setNames(factor(data[,var1]), data$species)
  y <- setNames(factor(data[,var2]), data$species)
  
  # ---- abor-single-levels ----
  if(any(length(levels(x)) == 1, 
         length(levels(y)) == 1)){
    return(data.frame(var1 = var1, var2 = var2, n = length(data$species), 
                      Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                      physig.var1 = NA, physig.var2 = NA,
                      cor = NA, phylocor = NA, 
                      lambda = NA, p = NA, ci = NA, 
                      error = paste("single state in", 
                                    if(length(levels(x)) == 1){var1},
                                    if(length(levels(y)) == 1){var2})))
  }
  
  # ---- get-phylosig ----
  physig.var1 <- eval(parse(
    text = paste(
      "phylo.d(data = data, phy = tree, names.col = species, binvar =", var1,")")))$DEstimate
  physig.var2 <- eval(parse(
    text = paste(
      "phylo.d(data = data, phy = tree, names.col = species, binvar =", var2,")")))$DEstimate
  
  # ---- get-model ----
  result.pgl <- try(fitPagel(tree = tree, x = x, y = y, method="fitMk"), silent = T)
  
  # ---- return-model ----
  if(result == "model"){
    return(results.pgl)
  }
  
  # ---- return-row ----
  if(class(result.pgl) == "try-error"){
    phylocor2 <- NA
    lambda <- NA
    p <- NA
    ci <- NA
    error <- gsub(pattern = ".*\n  ", "", geterrmessage())
  }else{
    phylocor2 <- result.pgl$lik.ratio
    lambda <- NA
    p <- result.pgl$P
    ci <- NA
    error <- NA
  }
  return(data.frame(var1 = var1, var2 = var2, n = length(data$species),
                    Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                    physig.var1, physig.var2,
                    cor = cor, phylocor = phylocor2, 
                    lambda = lambda, p = p, ci = ci, 
                    error = error))
}

fitCatDat <- function(data, tree, var1, var2, TD, result = "row") {
  
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$synonyms))
  tree$tip.label <- data$species[match(tree$tip.label, data$synonyms)]  
  
  # Std. correlation
  cor <- NA
  
  
  k <- length(levels(data[,var1]))
  I <- diag( k-1 )
  J <- matrix( rep(1, (k-1)^2), c(k-1, k-1))
  IJ <- (1/3) * (diag(k-1) + J)
  
  
  prior.phyl = list(R = list(V = IJ, fix = 1),G = list( G1 = list(V = IJ,
                                                                  n = k-1 ) ) )
  Ainv<-inverseA(tree)$Ainv
  
  results.mcmc <- try(MCMCglmm(as.formula(paste(var1, "~ trait + trait:", var2, "-1")), 
                               random = ~us(trait):species,
                               rcov = ~us(trait):units,
                               data = data,
                               ginverse = list(species=Ainv),
                               family = "categorical",
                               prior=prior.phyl,
                               thin   = 50,
                               burnin = 3000,
                               nitt   = 100000),
                      silent = F)
  if(result == "model"){
    return(results.mcmc)
  }
  
}

getLevels <- function(metadata, var) {
  meta <- list(levels = metadata[metadata$code == var, "scores"], 
               labels = metadata[metadata$code == var, "levels"]) %>%
    sapply(FUN = function(x){strsplit(x, ";")}, simplify = T) %>% 
    data.frame(stringsAsFactors = F)
  
  return(meta)
}

getTD <- function(data, tree, all.vars = T) {
  
  ptree <- drop.tip(tree, setdiff(tree$tip.label, data$synonyms))
  ptree$tip.label <- data$species[match(ptree$tip.label, data$synonyms)] 
  dist <- cophenetic(ptree)
  dspp <- rownames(dist)
  if(all.vars == T){
    comm.spp <- data$species[complete.cases(data)]
    comm <- setNames(data.frame(t(dspp %in% comm.spp)), dspp)}else{
      comm <- NULL
      vars <- names(data)[!names(data) %in% c("species", "synonyms")]
      for(var in vars){
        print(paste("calculating TD: ",which(vars == var),
                    "/", length(vars), "___", var,sep =""))
        comm.spp <- data$species[complete.cases(data[, var])]
        comm.r <- setNames(data.frame(t(dspp %in% comm.spp)), dspp)
        comm <- rbind(comm, comm.r)
      }}
  TD <-  taxondive(comm, dist)
  return(TD)
}


extractTD <- function(TD, vars) {
  sum.sts <- c("sd(Delta+)", "z(Delta+)",  "Pr(>|z|)")
  nvar <- length(vars)
  if(nvar == 1){
    sum.td <- data.frame(t(summary(td)[1,sum.sts]))
    TD <- TD %>% unlist() %>% t() %>% data.frame()
    TDdf <- cbind(var = vars, TD, sum.td)  
  }else{
    sum.td <- data.frame(summary(td)[1:nvar,sum.sts])
    TD <- do.call(cbind, TD)
    tdf <- cbind(var = vars, TD, sum.td)
  }
  return(tdf)
}


m1DataPrep <- function(g, phylo.match) {
  
  dat <- g$data
  datType <- unique(na.omit(g$meta_types[names(dat)]))
  char.err <-  sapply(dat, FUN = is.numeric)    
  
  if(!all(char.err)){
    print(paste("variables:", paste(names(dat)[!char.err], collapse = ", "),
                "contain non numeric values"))
    stop()}
  
  if("Int" %in% datType){
    int.vars <- intersect(names(dat), 
                          names(g$meta_types)[g$meta_types == "Int"]) 
    int.dat <- sapply(dat[,names(dat) %in% int.vars], FUN = is.integer)
    if(all(int.dat)){}else{
      for(int.var in int.vars){
        dat[,int.var] <-  as.integer(dat[,int.var])
      }
      print(paste("integer variables:", paste(int.vars, collapse = ", "),
                  "contain non integers. Rounded"))
    }
  }
  
  #Remove duplicate species matching to the same species on the tree
  dat <- dat[g$spp.list$species %in% phylo.match$species[
    phylo.match$data.status != "duplicate"],]
  # add synonym column to data 
  g$spp.list$synonyms <- phylo.match$synonyms[match(g$spp.list$species, 
                                                    phylo.match$species)]
  
  return(g)
  
  
}
## Create grid of unique variable combinations, calculate data availability for
## each and sort. min.n needs to be set in working environment.
varGridGen <- function(dat, ...) {
  
  var.grid <- calcTraitPairN(dat)
  var.grid <- var.grid[var.grid$n >= min.n,]
  var.grid <- var.grid[order(var.grid$n, decreasing = T),] 
  
  return(var.grid)
  
}

validateLog <- function(var, log.vars, data) {
  var %in% log.vars & all(na.omit(data[,var]) > 0)
}

validateInt <- function(var, data) {
  var %in% log.vars & is.integer(all(na.omit(data[,var])))
}

logData <- function(data, var) {
  data[,var] <- log(data[,var])
  names(data)[names(data) == var] <- paste("log", var, sep = "_")
  return(data)
}


getSppRanks <- function(ms_vars, data, load = T, input.folder, an.ID) {
  
  if(load & file.exists(paste(input.folder, "r data/", 
                              an.ID, "species_ranks.Rdata",
                              sep = ""))){
    load(paste(input.folder, "r data/", an.ID, "species_ranks.Rdata",
               sep = ""))
  }else{
    spp.vg <- expand.grid(ms_vars, ms_vars, stringsAsFactors = F)
    spp.vg <- spp.vg[apply(spp.vg, 1, FUN = function(x){length(unique(x)) > 1}),]
    
    spp.cs <- data.frame(species = data$species)
    
    for(i in 1:nrow(spp.vg)){
      dat.cs <- data[, unlist(spp.vg[i,])]
      spp.cs <- cbind(spp.cs, data.frame(complete.cases(dat.cs))) 
    }
    spp.rank <- data.frame(spp.cs$species, row.sums = rowSums(spp.cs[,-1]))
    spp.rank <- spp.rank[order(spp.rank$row.sums, decreasing = T),]
    save(spp.rank, spp.cs, spp.vg, 
         file = paste(input.folder, "r data/", an.ID, "species_ranks.Rdata",
                      sep = ""))
  }
  return(spp.rank)
}


logG <- function(g) {
  for(var in names(g$data)[names(g$data) %in% g$log.vars]){
    if(validateLog(var, g$log.vars, data = g$data)){
      g$data <- logData(g$data, var)
    for(elem in c("mgm_types", "meta_types", "lev")){
    names(g[[elem]])[names(g[[elem]]) == var] <- paste("log_", var, sep = "")}
    g$meta$code[g$meta$code == var] <- paste("log_", var, sep = "")
    }
  }
  return(g)
}



testImpute <- function(data, test.params, ms_vars, metadata, ...) {
  
  error <- NULL
  
  for(i in 1:nrow(test.params)){
    v.p <- test.params$vars[i]
    s.p <- test.params$spp[i]
    
    t0 <- Sys.time()
    print(paste(i, "/", nrow(test.params),"|_____", 
                paste(test.params[i,], collapse = "-"), sep = ""))
    
    vars <- ms_vars[1:(as.integer(length(ms_vars)*v.p))]
    spps <- spp.ranks$spp[1:(as.integer(nrow(spp.ranks)*s.p))]
    spp <- data$species %in% spps
    
    g <- list(data = data[spp, vars], 
              meta = metadata, 
              mgm_types = mgm_types, meta_types = meta_types, 
              log.vars = log.vars, 
              spp.list = data.frame(species = data[spp, c("species")]),
              v.p = v.p, s.p = s.p)
    
    source(paste(script.folder, "mgm_dataprep.R", sep = ""))
    
    #' impute-mgm
    source(paste(script.folder, "impute_mgm.R", sep = ""))
    impt <- g[c("imp_data","OOBerror")]
    impt.out <- c(impt, list(vars = vars, spps = spps, v.p = v.p, s.p = s.p))
    save(impt.out, file = paste(output.folder, "data/imputed_data/", an.ID, 
                                "-v", test.params$vars[i], 
                                "-s", test.params$spp[i],
                                ".Rdata", sep = ""))
    error <- rbind(error, c(g$OOBerror, v.p = v.p, s.p = s.p))
    print(Sys.time() - t0)
  }
  error <- data.frame(error)
  error.s <- error[order(error$NRMSE),]
  write.csv(error.s, paste(input.data, "r data/imp_err_", an.ID, ".Rdata",
                           sep = ""))
  return(error.s)
}



