#' calculate number of variable pairwise data-point availability 
#'
#' @param data a wide data_frame containing variables in columns and species in rows. Function ignores columns named "species" and "synonyms". 
#' @details 
#' 
#' @return a data.frame of the unique variable pairwise combinations and counts (n) of data availability. The dataframe is sorted according to $n highest to lowest
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

#' Title
#'
#' @param pair 
#' @param data 
#' @param tree 
#' @param log.vars 
#' @param pair.type 
#' @param result 
#' @param mgm_types 
#' @param alco 
#' @param overwrite 
#'
#' @return
#' @export
#'
#' @examples
pglsPhyloCor <- function(pair, data, tree, log.vars = NULL, pair.type, 
                         result = "row", mgm_types = mgm_types, alco = F, overwrite = F){
  pair <- as.character(pair)
  error <- NULL
  
  
  if(file.exists(paste0(script.folder, "data/phy_models/", an.ID, "_", 
                        if(alco){"alco_"}else{}, 
                        pair[1], "-", pair[2],".rda")) & !overwrite){
    load(file = paste0(script.folder, "data/phy_models/", an.ID, "_", 
                       if(alco){"alco_"}else{}, pair[1], "-", pair[2],".rda"))
    if(pair.type != "cc"){
      if(result == "row"){
        if(exists("data_row", envir = environment(), inherits = F)){
          if(!"alco_id" %in% names(data_row)){
            if(any(pair %in% alco.vars)){
              if(alco){
                alco_id <- "alco"
              }else{alco_id <- "non_alco"}
            }else{alco_id <- NA}
            data_row$alco_id <- alco_id
            }
          cat("returning loaded data_row", "\n")
          return(data_row)
        }}
      
      if(result == "model"){
        if(exists("result.pgls", envir = environment(), inherits = F)){
          cat("returning loaded result.pgls")
          return(result.pgls) 
        }
      }
    }else{
      cat("returning loaded tau. Use get_cc_row to extract data into row")
      return(tau)
    }
  }
  

  
  # ---- set-up-vars ----
  data <- data[, c("species", pair)]
  
  # ---- get-TD ----
  if(alco & pair.type != "cc"){
    TD <- NA
    class(TD) <- "ignore"
  }else{  
    TD <- getTD(data, tree)
  }
  
  # ---- sub-complete-cases ----
  data <- data[complete.cases(data),]
  
  # ---- aggregate-cats ----  
  for(var in pair){
    if(mgm_types[var] == "c"){
      
      data[,var] <- aggVar(data[,var], min.cat)
      data <- data[complete.cases(data),]
      data[,var] <- as.factor(data[,var])
      if(1 >= length(levels(data[,var]))){
        error <- list(pair = pair,
                    error = paste0("ERROR: singular category: '",  levels(data[,var]),
                                   "' in variable ", var))}
    }
  }
  
  # ---- fit-phylocor ---- 
  if(pair.type %in% c("nc", "nn")) {
    row <- fitNumDat(data, tree, pair, TD, result, mgm_types = mgm_types, 
                     error, pair.type, alco = alco)
  }
  if(pair.type == "cc"){
    row <-  fitGKtau(data, tree, pair, TD, result = "row", error)
  }
  if(pair.type == "bb") {
    row <-  fitBinDat(data, tree, pair, TD, result, error) 
  }
  if(pair.type == "ccMCMC") {
    meta <- list(getLevels(metadata, pair[1]), getLevels(metadata, pair[2]))
    # factorise data
    data[,c(pair)] <- mapply(FUN = function(x, meta){factor(as.character(x),
                                                                  levels = meta$levels,
                                                                  labels = meta$labels)},
                                   x = data[,pair], meta = meta) %>%
      as.data.frame(stringsAsFactors = T)
    
    row <-  fitCatDat(data, tree, pair, TD, result) }
  

  
  # ---- return ----
  return(row)
  
}

fitNumDat <- function(data, tree, pair, TD, result = "row", mgm_types = mgm_types, 
                      error, pair.type, alco) {
  
  if(any(pair %in% alco.vars)){
    if(alco){
      alco_id <- "alco"
    }else{alco_id <- "non_alco"}
  }else{alco_id <- NA}
  
  if(class(TD) == "ignore"){
    Dplus <- NA 
    sd.Dplus <- NA 
    EDplus <- NA
  }else{
    Dplus <- TD$Dplus 
    sd.Dplus <- TD$sd.Dplus 
    EDplus <- TD$EDplus
  }
    
  # ---- configure-tree ----
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))

  # ---- get-phylosig ----
  if(!alco){
    physig.var1 <- phylosig(tree, setNames(data[,pair[1]], data$species),
                            method = "lambda")$lambda
    }else{
      physig.var1 <- NA
    }
  
  if(mgm_types[pair[2]] == "g"){
    if(!alco){
      physig.var2 <- phylosig(tree, setNames(data[,pair[2]], data$species),
                              method = "lambda")$lambda
    }else{
      physig.var2 <- NA  
    }
  }else{physig.var2 <- NA}
  
  if(!is.null(error)){
    class(error) <- "try-error"
    result.pgls <- error
    if(!alco){
      data_row <- data.frame(var1 = pair[1], var2 = pair[2], n = length(data$species),
                             Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                             physig.var1, physig.var2,
                             cor = NA, phylocor = NA, 
                             lambda = NA, p = NA, l.ci.l = NA, l.ci.u = NA,
                             error = error$error, pair.type = pair.type, alco_id = alco_id)
    }else{
      data_row <- data.frame(var1 = pair[1], var2 = pair[2], n = length(data$species),
                              Dplus = Dplus, sd.Dplus = sd.Dplus, EDplus = EDplus,
                              physig.var1, physig.var2,
                              cor = NA, phylocor = NA, 
                              lambda = NA, p = NA, l.ci.l = NA, l.ci.u = NA,
                              error = error$error, pair.type = pair.type, 
                              alco_id = alco_id)
    }
    save(result.pgls, data_row, file = paste0(script.folder, "data/phy_models/", 
                                              an.ID, "_", if(alco){"alco_"}else{}, 
                                              pair[1], "-", pair[2],".rda"))
    if(result == "model"){
      return(result.pgls)}else{
        return(data_row)
      }
  }
  
  # ---- get-cor ----
  if(mgm_types[pair[2]] == "g"){
    if(!alco){
      cor <- cor(data[,pair[1]], data[,pair[2]])
    }else{
      cor <- NA 
    }}else{
    cor <- NA}
  
  ## ---- get-models ----
  if(exists("result.pgls", envir = parent.frame(), inherits = F)){
    result.pgls <- get("result.pgls", envir = parent.frame())
  }else{
    cd <- comparative.data(phy = tree, data = data, names.col = "species", vcv=F)
    result.pgls <- try(pgls(as.formula(paste(pair[1], "~", pair[2], sep = "")),
                            data = cd, lambda="ML"))
  }
  # ---- return-model ----
  if(result == "model"){
    return(result.pgls)
  }
  
  # ---- return-row ---- 
  if(class(result.pgls) == "try-error"){
    phylocor2 <- NA
    lambda <- NA
    p <- NA
    ci.l <- NA
    ci.u <- NA
    error <- gsub(pattern = ".*\n  ", "", geterrmessage())
  }else{
    if(mgm_types[pair[2]] == "g"){
      t <- summary(result.pgls)$coefficients[pair[2],"t value"]
      df <- as.vector(summary(result.pgls)$fstatistic["dendf"])
      phylocor2 <- sqrt((t*t)/((t*t)+df))*sign(summary(result.pgls)$coefficients[pair[2],"Estimate"])
      p <- coef(summary(result.pgls))[2,4]
    }else{
      phylocor2 <- sqrt(summary(result.pgls)$r.squared)
      p <- p_value(result.pgls)
    }
    lambda <- result.pgls$param["lambda"]
    ci.l <- result.pgls$param.CI$lambda$ci.val[1]
    ci.u <- result.pgls$param.CI$lambda$ci.val[2]
    error <- NA
  }
  
  # ---- data-row ----
  data_row <- data.frame(var1 = pair[1], var2 = pair[2], n = length(data$species),
                         Dplus = Dplus, sd.Dplus = sd.Dplus, EDplus = EDplus,
                         physig.var1, physig.var2,
                         cor = cor, phylocor = phylocor2, 
                         lambda = lambda, p = p, l.ci.l = ci.l, l.ci.u = ci.u,
                         error = error, pair.type = pair.type, alco_id = alco_id)
  # ---- save ----
  save(result.pgls, data_row, file = paste0(script.folder, "data/phy_models/", an.ID, "_", 
                                  if(alco){"alco_"}else{}, pair[1], "-", pair[2],".rda"))
  return(data_row)

}

fitBinDat <- function(data, tree, pair, TD, result = "row") {

  # ---- configure-tree ----
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
  tree$tip.label <- data$species[match(tree$tip.label, data$species)]  
  
  # ---- get-cor ----
  cor <- NA
  
  # ---- get-cor ----
  data <- data[match(tree$tip.label, data$species),]
  x <- setNames(factor(data[,pair[1]]), data$species)
  y <- setNames(factor(data[,pair[2]]), data$species)
  
  # ---- abor-single-levels ----
  if(any(length(levels(x)) == 1, 
         length(levels(y)) == 1)){
    return(data.frame(var1 = pair[1], var2 = pair[2], n = length(data$species), 
                      Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                      physig.var1 = NA, physig.var2 = NA,
                      cor = NA, phylocor = NA, 
                      lambda = NA, p = NA, l.ci.l = NA, l.ci.u = NA,
                      error = paste("single state in", 
                                    if(length(levels(x)) == 1){pair[1]},
                                    if(length(levels(y)) == 1){pair[2]})),
           pair.type = "bin")
  }
  
  # ---- get-phylosig ----
  physig.var1 <- eval(parse(
    text = paste(
      "phylo.d(data = data, phy = tree, names.col = species, binvar =", pair[1],")")))$DEstimate
  physig.var2 <- eval(parse(
    text = paste(
      "phylo.d(data = data, phy = tree, names.col = species, binvar =", pair[2],")")))$DEstimate
  
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
  return(data.frame(var1 = pair[1], var2 = pair[2], n = length(data$species),
                    Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                    physig.var1, physig.var2,
                    cor = cor, phylocor = phylocor2, 
                    lambda = lambda, p = p, ci = ci, 
                    error = error, pair.types = "bb"))
}

fitCatDat <- function(data, tree, pair, TD, result = "row") {
  
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
  tree$tip.label <- data$species[match(tree$tip.label, data$species)]  
  
  # Std. correlation
  cor <- NA
  
  
  k <- length(levels(data[,pair[1]]))
  I <- diag( k-1 )
  J <- matrix( rep(1, (k-1)^2), c(k-1, k-1))
  IJ <- (1/3) * (diag(k-1) + J)
  
  
  prior.phyl = list(R = list(V = IJ, fix = 1),G = list( G1 = list(V = IJ,
                                                                  n = k-1 ) ) )
  Ainv<-inverseA(tree)$Ainv
  
  results.mcmc <- try(MCMCglmm(as.formula(paste(pair[1], "~ trait + trait:", pair[2], "-1")), 
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

fitGKtau <- function(data, tree, pair, TD, result = "row", error) {
  
  if(!is.null(error)){
    error$TD <- TD
    return(error)}
 
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
  tau <- get_simGKtau(pair, data, tree, nsim = 500)
  tau <- list(pair = pair, phy.GKtau = tau, GKtau = GKtau(data[, pair[1]], data[, pair[2]]),
              TD = TD, error = NA)
  save(tau, file = paste0(script.folder, "data/phy_models/", an.ID, "_", 
                          pair[1], "-", pair[2],".rda"))
   return(tau)
} 

get_cc.row <- function(x, cc_vg){
  
  if(length(x) == 2){
    row <- data.frame(var1 = x$pair[1], var2 = x$pair[2], 
                      n = cc_vg$n[apply(cc_vg[,1:2], 1, 
                                        FUN = function(y){identical(unname(y), x$pair)})], 
                      Dplus = NA, sd.Dplus = NA, EDplus = NA,
                      physig.var1 = NA, physig.var2 = NA, cor = NA,  
                      phylocor = NA, lambda = NA, p = NA, l.ci.l = NA, 
                      l.ci.u = NA, 
                      error = x$error,
                      pair.type = "cc", alco_id = NA)
    return(row)
  }
  
  if(length(x) == 3 & !is.null(x$TD)){
    row <- data.frame(var1 = x$pair[1], var2 = x$pair[2], 
                      n = cc_vg$n[apply(cc_vg[,1:2], 1, 
                                        FUN = function(y){identical(unname(y), x$pair)})], 
                      Dplus = x$TD$Dplus, sd.Dplus = x$TD$sd.Dplus, EDplus = x$TD$EDplus,
                      physig.var1 = NA, physig.var2 = NA, cor = NA,  
                      phylocor = NA, lambda = NA, p = NA, l.ci.l = NA, 
                      l.ci.u = NA, 
                      error = x$error,
                      pair.type = "cc", alco_id = NA)
    return(row)
  }
  
  require(boot)
  
  max.id <- which(x$GKtau[,5:6] == max(x$GKtau[,5:6]))[1]
  min.id <-setdiff(1:2, max.id)
  obs <-  unlist(x$GKtau[,5:6][max.id], use.names = F)
  simbase <- unlist(x$phy.GKtau)
  to_NA <- as.logical(is.infinite(simbase) + is.nan(simbase))
  simbase <- na.omit(replace(simbase,to_NA, NA))
  
  ci <- boot.ci(boot(simbase,function(x,i) median(x[i]), R=1000), type = "basic")
  base <- ci$t0
  l.ci.l <- ci$basic[4]
  l.ci.u <- ci$basic[5]
  
  if(base > obs){
    phylocor <- 0
  }else{
    phylocor <- (obs-base)/(1-base)
    }
 
  n <- length(simbase)
  r <- sum(simbase > obs)
  p <- (r+1)/(n+1)
  
  row <- data.frame(var1 = x$pair[max.id], var2 = x$pair[min.id], 
                    n = cc_vg$n[apply(cc_vg[,1:2], 1, 
                                      FUN = function(y){identical(unname(y), x$pair)})], 
                    Dplus = x$TD$Dplus, sd.Dplus = x$TD$sd.Dplus, EDplus = x$TD$EDplus,
                    physig.var1 = NA, physig.var2 = NA, cor = max(x$GKtau[,5:6]),  
                    phylocor = phylocor, lambda = base, p = p, l.ci.l, l.ci.u, error = NA,
                    pair.type = "cc", alco_id = NA)
  return(row)
}


get_traitQ <- function(trait, tree) {
  trait.er<-rerootingMethod(tree, trait, model="ER")
  trait.sym<-rerootingMethod(tree, trait, model="SYM")
  if(pchisq(-2*(trait.er$loglik-trait.sym$loglik), df=3-1, lower.tail=FALSE) > 0.05){
    Q <- trait.er$Q
    attr(Q, which = "model") <- "ER"
  }else{
    Q <- trait.sym$Q
    attr(Q, which = "model") <- "SYM"
  }
  return(Q)
}

get_traitSim <- function(Q, tree, nsim = 500) {
  sims <- as.data.frame(sim.char(tree, par = Q, nsim = nsim, model="discrete"))
  return(sims)
}

get_simGKtau <- function(pair, data, tree, nsim = 500) {
  
  Q1 <- get_traitQ(setNames(data[,pair[1]], data$species), tree)
  sims1 <- get_traitSim(Q = Q1, tree, nsim = nsim)
  Q2 <- get_traitQ(setNames(data[,pair[2]], data$species), tree)
  sims2 <- get_traitSim(Q = Q2, tree, nsim = nsim)
  
  tau <- mapply(FUN = function(x, y){
    tau <- GKtau(x,y)
    tau[, c("tauxy", "tauyx")]}, x = sims1, y = sims2, SIMPLIFY = F) %>% do.call(what = rbind)
  attr(tau, "model") <- c(attr(Q1, "model"), attr(Q2, "model"))
  
  return(tau)
}

getLevels <- function(metadata, var) {
  meta <- list(levels = metadata[metadata$code == var, "scores"], 
               labels = metadata[metadata$code == var, "levels"]) %>%
    sapply(FUN = function(x){strsplit(x, ";")}, simplify = T) %>% 
    data.frame(stringsAsFactors = F)
  
  return(meta)
}

getTD <- function(data, tree, all.vars = T) {
  
  ptree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
  ptree$tip.label <- data$species[match(ptree$tip.label, data$species)] 
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


validate_log.vars <- function(log.vars, data) {
  valid <- NULL
  for(i in 1:length(log.vars)){
    var <- log.vars[i]
    valid <- c(valid, all(na.omit(data[,var]) > 0))
  }
  names(valid) <- log.vars
  valid
}



validateInt <- function(var, data) {
  var %in% log.vars & is.integer(all(na.omit(data[,var])))
}

corr_zero.logs <- function(data, replace, replace.with = NA) {
  for(i in 1:length(replace)){
    var <- replace[i]
    data[which(data[,var] == 0), var] <- replace.with
  }
  data
}

logData <- function(data, var) {
  data[,var] <- log(data[,var])
  names(data)[names(data) == var] <- paste("log", var, sep = "_")
  return(data)
}


getSppRanks <- function(ms_vars, data, load = T, input.folder, an.ID) {
  
  if(load & file.exists(paste(input.folder, "taxo/", 
                              an.ID, "species_ranks.Rdata",
                              sep = ""))){
    load(paste(input.folder, "taxo/", an.ID, "species_ranks.Rdata",
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
         file = paste(input.folder, "taxo/", an.ID, "species_ranks.Rdata",
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


#' Title
#'
#' @param vg 
#' @param mgm_types 
#'
#' @return
#' @export
#'
#' @examples
get_vg_dt <- function(vg, mgm_types) {
  vg_dt <- data.frame(Var1 = mgm_types[vg$Var1], Var2 = mgm_types[vg$Var2])
  return(vg_dt)
}


#' Title
#'
#' @param vg 
#' @param mgm_types 
#'
#' @return
#' @export
#'
#' @examples
orderby_vg_dt <- function(vg, mgm_types){
  
  vg_dt <- get_vg_dt(vg, mgm_types)
  cv <- as.data.frame(t(apply(vg_dt, 1, function(x){order(x, decreasing = T)})))
  t.vg <- vg
  t.vg$Var1 <- vg[cbind(1:nrow(vg),cv$V1)]
  t.vg$Var2 <- vg[cbind(1:nrow(vg),cv$V2)]
  vg <- t.vg
  return(vg)
}




#' sex_data_merge
#'
#' @param data 
#' @param proc.pars 
#'
#' @return
#' @export
#'
#' @examples
sex_data_merge <- function(data, proc.pars) {
  # supporting functions
  get_ssd <- function(male, female){
    log(ds$data[,male]/ds$data[,female])
  }
  merge_means <- function(name, primary, mean, male, female) {
    sex.means <- rowMeans(ds$data[, c(female, male)]) %>%
      setNames(ds$data$species)
    if(length(mean) != 0){
      means <- ds$data[,mean] %>%
        setNames(ds$data$species)
    }
    
    # set primary & supplement with secondary source of data (if applicable)
    if(primary == "mean"){ # primary = mean
      comp <- means
      add <- names(sex.means) %in% names(na.omit(sex.means))[!names(na.omit(sex.means)) %in% names(na.omit(means))]
      comp[add] <- sex.means[add]
    }else{ # primary = sex.means
      comp <- sex.means 
      if(length(mean) != 0){
        add <- names(means) %in% names(na.omit(means))[!names(na.omit(means)) %in% names(na.omit(sex.means))]
        comp[add] <- means[add]
      }
    } 
    # remove source vars. add comp data  
    ds$data <- ds$data[,!names(ds$data) %in% c(mean, female, male)]
    ds$data[,name] <- comp
  }
  process_nest.row <- function(x){
    primary = x$primary
    name = x$name
    table <- x$data[[1]]
    mean <- table[table$sex == "mean", "code"] %>% .[["code"]]
    male <-  table[table$sex == "male", "code"] %>% .[["code"]]
    female <- table[table$sex == "female", "code"] %>% .[["code"]]
    
    cat(paste0("processing comp_var: ", name, "\n"))
    # append ssd [log(male/female)] to data
    ds$data[,paste0(x$name, ".ssd")] <- get_ssd(male, female)
    
    # merge sex data and means
    merge_means(name, primary, mean, male, female)
  }
  
  #create data environment
  ds <- new.env(parent = emptyenv())
  assign("data", data, envir = ds)
  
  nest.pars <- proc.pars %>% 
    filter(!is.na(name)) %>% 
    group_by(name, primary) %>% 
    nest()
  
  for(i in 1:nrow(nest.pars)){
    process_nest.row(nest.pars[i,])
  }
  
  ds$data <- select(ds$data, species:order, everything()) 
  return(ds$data)
}


get_an_meta <- function(metadata, data, proc.pars){
  ammend_comb.var_descr <- function(comb.vars, an_meta) {
    var_ids <- an_meta$code %in% comb.vars
    an_meta$descr[var_ids] <- paste(an_meta$descr[var_ids], "(mean + sex combined)")
    return(an_meta)
  }
  get_new_comb_meta_rows <- function(new_comb.var, metadata, proc.pars) {
    get_new_comb_meta_row <- function(new_comb.var, metadata, proc.pars) {
      row_var <- proc.pars$code[proc.pars$name %in% new_comb.var][1]
      meta_row <- metadata[metadata$code == row_var,]
      meta_row$code <- new_comb.var
      meta_row$descr <- gsub("male|female", "sex combined", meta_row$descr)
      return(meta_row)
    }
    do.call(rbind, lapply(new_comb.vars, FUN = get_new_comb_meta_row, 
                          metadata = metadata, proc.pars = proc.pars))
  }
  get_ssd_meta_rows <- function(ssd.vars, an_meta, proc.pars) {
    add <- an_meta[match(gsub(".ssd", "", ssd.vars), an_meta$code),]
    add$code <- ssd.vars
    add$descr <- paste(add$descr, "ssd [log(m/f)]")
    add$log <- F
    add$units <- NA
    return(add)
  }
  
  an_meta <- metadata[metadata$code %in% names(data)[
    names(data) != "species"],]
  
  add_vars <- names(data)[!names(data) %in% an_meta$code]
  new_comb.vars <- unique(proc.pars$name[proc.pars$name %in% add_vars])
  comb.vars <- unique(na.omit(proc.pars$name))[
    !unique(na.omit(proc.pars$name)) %in% new_comb.vars]
  ssd.vars <- add_vars[grep(".ssd", add_vars)]
  
  an_meta <- rbind(ammend_comb.var_descr(comb.vars, an_meta),
                   get_new_comb_meta_rows(new_comb.vars, metadata, proc.pars))
  an_meta <- rbind(an_meta, get_ssd_meta_rows(ssd.vars, an_meta, proc.pars)) %>% arrange(code)
  
  return(an_meta)
}

alco_correct_data <- function(data, alco.vars, tree, mgm_types, pair.type){
  alco.complete <- setNames(alco.vars %in% names(data), alco.vars)
  alco.vars <- alco.vars[alco.complete]
  
  for(alco.var in alco.vars){
    cat("allometric correcting:", alco.var, "(", which(alco.vars == alco.var), "of", length(alco.vars), ") -")
    result <- pglsPhyloCor(pair = c(alco.var,"body.mass_log"), data = data, tree = tree, log.vars = log.vars, 
                           pair.type = pair.type, result = "model", mgm_types = mgm_types, alco = T)
    if(class(result)[1] == "try-error"){
      alco.complete[alco.var] <- F 
      cat("not successful", "\n")
    }else{
      data[, alco.var] <- NA
      data[match(rownames(result$residuals), data$species), alco.var] <- result$phyres
      cat("successful", "\n")
    }
  }
  return(list(data = data, alco.complete = alco.complete))
}



#' run rnetcarto analysis
#'
#' The function runs a network analysis on the phylocor results and saves and returns
#' required outputs and statistics
#' @param res results data.frame from the PhyloNetworker analysis
#' @param edge_det character string. The edge determination method. One of "p-value" 
#' or "phylocor"
#' @param alco logical. Whether allometric correction has been applied. Ideally 
#' specified through param file.
#' @param save logical. Whether results should be saved to disk at path [`output.folder`/"data/networks/"]
#'
#' @return a data.frame containing modified rnetcarto output, appended with additional metadata. If save = T, 
#' an additional data.frame containing network level statistics is also written to disk
#' @export
#'
#' @examples
run_netcarto <- function(res, edge_det, alco, save = F) {
  # real modularity function
  real_modularity <- function(net.d, net_tbl, metric = 2) {
    q_r <- cbind(net_tbl$module[match(net.d$var1, net_tbl$name)], 
                 net_tbl$module[match(net.d$var2, net_tbl$name)]) %>% 
      apply(1, FUN = function(x){length(unique(x)) == 1}) %>% 
      sum()/nrow(net.d)
    if(metric == 1){return(q_r)}else{return(2 * q_r -1)}
  }
  
  # check edge_det argument
  if(!edge_det %in% c("p-value", "phylocor")){
    stop('edge_det must be one of c("p-value", "phylocor")')
  }
  # ---- pn-sub-res ----
  if(edge_det == "p-value"){
    res_sub <- res[res[,"p"] < 0.05 & !is.na(res[,"p"]),]
    net.d <- data.frame(res_sub[, c("var1", "var2")],
                        edge_det = res_sub[,"p"])
  }else{
    res_sub <- res[abs(res$phylocor) > cutoff & !is.na(res$phylocor),]
    net.d <- data.frame(res_sub[, c("var1", "var2")],
                        edge_det = res_sub[,"phylocor"])
  }
  # ---- pn-network ----
  net.list <- as.list(net.d)
  net <- netcarto(web = net.list, seed = 1)
  net_tbl <- cbind(net[[1]], modularity = net[[2]])
  
  # ---- pn-net_stats ----
  nvar <- length(unique(unlist(res[,1:2])))
  tpairs <- choose(nvar,2) 
  errors <- sum(!is.na(res$error)) - length(grep("ERROR: singular category:", res$error))
  pairs <- nrow(res) - length(grep("ERROR: singular category:", res$error))
  edges.phy <- dim(res[abs(res$phylocor) > cutoff & !is.na(res$phylocor),])[1]
  edges.p <- dim(res[res[,"p"] < 0.05 & !is.na(res[,"p"]),])[1]
  nodes <- nrow(net_tbl)
  q_r2 <- real_modularity(net.d, net_tbl, metric = 2)
  net_stats <- data.frame(nvar, tpairs, errors, pairs, edges.phy, 
                          edges.p, nodes, q_r2)
  
  net_tbl <- cbind(net_tbl, an_meta[match(net_tbl$name, an_meta$code), c("descr", "cat", "type", "scores", "levels")])
  net_tbl <- cbind(net_tbl, TDdf[match(net_tbl$name, TDdf$var), c("EDplus", "Dplus", "sd.Dplus", "phylosig")])
  
  net_tbl$p.prop <- NA
  net_tbl$pairs_tested <- NA
  net_tbl$phycor_mean <- NA
  net_tbl$phycor_median <- NA
  net_tbl$phycor_max <- NA
  
  for(node in net_tbl$name){
    node.rows <- as.logical((res$var1 == node) + (res$var2 == node))
    ps <- na.omit(res[node.rows,"p"])
    cs <- abs(na.omit(res[node.rows,"phylocor"]))
    p.pairs <- length(ps)
    net_tbl$pairs_tested[net_tbl$name == node] <- p.pairs
    net_tbl$p.prop[net_tbl$name == node] <- sum(ps < 0.05)/p.pairs
    net_tbl$phycor_mean[net_tbl$name == node] <- mean(cs)
    net_tbl$phycor_median[net_tbl$name == node] <- median(cs)
    net_tbl$phycor_max[net_tbl$name == node] <- max(cs)
  }
  
  if(save){
    write.csv(net_tbl, 
              paste(output.folder, "data/networks/", an.ID,
                    "_net_mn", min.n, 
                    "_", edge_det,
                    if(alco){"alco_"},
                    if(log){"_log"},".csv", sep = ""),
              row.names = F)
    write.csv(net_stats, 
              paste(output.folder, "data/networks/", an.ID,
                    "_net_mn", min.n, 
                    "_", edge_det,
                    if(alco){"alco_"},
                    if(log){"_log"},"_net-stats.csv", sep = ""),
              row.names = F)
    save(net, file = paste(output.folder, "data/networks/", an.ID,
                           "_net_mn", min.n, 
                           "_", edge_det,
                           if(alco){"alco_"},
                           if(log){"_log"},".RData", sep = ""))
  }
  return(net_tbl)
}


#' Get variable grid row indices
#'
#' @param vars character vector. 
#' @param vg a `var_grid` data.frame of trait pairs
#'
#' @return If `vars` contains string "all", an index for each 
#' row of the vg data.frame is supplied. If a vector of variable names is supplied to `vars`,
#' a vector of row indices of vg trait pairs containing any of the variables supplied is returned.
#' @export
#'
#' @examples
get_vgr_id <- function(vars, vg) {
  if("all" %in% vars ){1:nrow(vg)}else{
    sort(c(which(vg$Var1 %in% vars), which(vg$Var2 %in% vars)))
  }
}

