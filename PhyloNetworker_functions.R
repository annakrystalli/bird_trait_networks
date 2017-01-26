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

pglsPhyloCor <- function(pair, data, tree, log.vars = NULL, pair.type, 
                         result = "row", mgm_types = mgm_types){
  pair <- as.character(pair)
  error <- NULL
  
  # ---- set-up-vars ----
  data <- data[, c("species", pair)]
  
  # ---- get-TD ----  
  TD <- getTD(data, tree)
  
  
  # ---- sub-complete-cases ----
  data <- data[complete.cases(data),]
  
  # ---- log-log.vars ----  
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
    if(mgm_types[var] == "g"){
      if(validateLog(var = var, log.vars, data)){
        data[, var] <- log(data[, var])
        colnames(data)[colnames(data) == var] <- paste("log", var, sep = "_")
        pair[pair == var] <- paste("log", var, sep = "_")
        names(mgm_types)[names(mgm_types) == var] <- paste("log", var, sep = "_")
      }
    }
  }
  
  # ---- fit-cor ---- 
  if(pair.type %in% c("nc", "nn")) {
    row <- fitNumDat(data, tree, pair, TD, result, mgm_types = mgm_types, 
                     error, pair.type)
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
                      error, pair.type) {
  
  # ---- configure-tree ----
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))

  # ---- get-phylosig ----
  physig.var1 <- phylosig(tree, setNames(data[,pair[1]], data$species),
                          method = "lambda")$lambda
  if(mgm_types[pair[2]] == "g"){
    physig.var2 <- phylosig(tree, setNames(data[,pair[2]], data$species),
                            method = "lambda")$lambda
  }else{
    physig.var2 <- NA  
  }
  
  if(!is.null(error)){
    return(data.frame(var1 = pair[1], var2 = pair[2], n = length(data$species),
                      Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                      physig.var1, physig.var2,
                      cor = NA, phylocor = NA, 
                      lambda = NA, p = NA, l.ci.l = NA, l.ci.u = NA,
                      error = error$error, pair.type = pair.type))
  }
  
  
  # ---- get-cor ----
  if(mgm_types[pair[2]] == "g"){
    cor <- cor(data[,pair[1]], data[,pair[2]])
  }else{
    cor <- NA 
  }
  
  ## ---- get-model ----
  cd <- comparative.data(phy = tree, data = data, names.col = "species", vcv=F)
  
  result.pgls <- try(pgls(as.formula(paste(pair[1], "~", pair[2], sep = "")),
                          data = cd, lambda="ML"))
  
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
  return(data.frame(var1 = pair[1], var2 = pair[2], n = length(data$species),
                    Dplus = TD$Dplus, sd.Dplus = TD$sd.Dplus, EDplus = TD$EDplus,
                    physig.var1, physig.var2,
                    cor = cor, phylocor = phylocor2, 
                    lambda = lambda, p = p, l.ci.l = ci.l, l.ci.u = ci.u,
                    error = error, pair.type = pair.type))
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
                      pair.type = "cc")
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
                      pair.type = "cc")
    return(row)
  }
  
  require(boot)
  
  max.id <- which(x$GKtau[,5:6] == max(x$GKtau[,5:6]))[1]
  min.id <-setdiff(1:2, max.id)
  obs <-  unlist(x$GKtau[,5:6][max.id], use.names = F)
  sim.base <- na.omit(replace(x$phy.GKtau[,max.id],
                              is.infinite(x$phy.GKtau[,max.id]),
                              NA))
  
  ci <- boot.ci(boot(sim.base,function(x,i) median(x[i]), R=1000), type = "basic")
  base <- ci$t0
  l.ci.l <- ci$basic[4]
  l.ci.u <- ci$basic[5]
  
  if(base > obs){phylocor <- 0}else{
    phylocor <- (obs-base)/(1-base)}
  
  row <- data.frame(var1 = x$pair[max.id], var2 = x$pair[min.id], 
                    n = cc_vg$n[apply(cc_vg[,1:2], 1, 
                                      FUN = function(y){identical(unname(y), x$pair)})], 
                    Dplus = x$TD$Dplus, sd.Dplus = x$TD$sd.Dplus, EDplus = x$TD$EDplus,
                    physig.var1 = NA, physig.var2 = NA, cor = max(x$GKtau[,5:6]),  
                    phylocor = phylocor, lambda = base, p = NA, l.ci.l, l.ci.u, error = NA,
                    pair.type = "cc")
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

validateInt <- function(var, data) {
  var %in% log.vars & is.integer(all(na.omit(data[,var])))
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