library(caper)
library(geiger)
library(PHYLOGR)


# produce named vector of variable data to use for phyloCor analysis
getNamedVector <- function(var, data){
  x <- as.vector(as.numeric(data[,var]))
  names(x) <- data$species
  return(x)
}

# Produce new phylo.matrix for subset of species
subsetPhylomat <- function(spp, phylomat, match.dat = NULL){
  
  nsps <- length(spp) 
  vmat <- matrix(NA, nrow = nsps, ncol = nsps, dimnames = list(spp, spp))
  
  m.id <- cbind(rep(spp, times = nsps), rep(spp, each = nsps))
  
  if(is.null(match.dat)){
    if(all(spp %in% unlist(dimnames(phylomat)))){p.id <- m.id}else{
      stop("data species names do not match phylogeny tip names. correct or provide match.dat data.frame")}
  }else{
    p.id <- cbind(rep(match.dat$synonyms[match(spp, match.dat$species)], times = nsps), 
                  rep(match.dat$synonyms[match(spp, match.dat$species)], each = nsps))}
  
  vmat[m.id] <- phylomat[p.id]
  
  return(vmat)
}

phylo.mean <- function(x, phylomat){
  mean <-colSums(phylomat%*%x)/sum(phylomat)}

phylo.var <- function(x, mean, phylomat, nsps){
  var.x <-t(x-mean) %*% phylomat%*%(x-mean)/(nsps-1)
}

getPhyloCor <- function(x, y, phylomat, nsps){
  
  if(any(names(x) != names(y))){stop("vector species names mismatch")}
  if(dim(phylomat)[1] != dim(phylomat)[2]){stop("phylomat not square")}
  if(any(dimnames(phylomat)[[1]] != dimnames(phylomat)[[2]])){stop("phylomat dimnames mismatch")}
  if(any(names(x) != dimnames(phylomat)[[1]], names(x) != dimnames(phylomat)[[2]])){stop("x and phylomat name mismatch")}
  if(any(names(y) != dimnames(phylomat)[[1]], names(y) != dimnames(phylomat)[[2]])){stop("y and phylomat name mismatch")}
  
  
  mean.x <- phylo.mean(x, phylomat)
  mean.y <- phylo.mean(y, phylomat)
  
  var.x <- phylo.var(x, mean.x, phylomat, nsps = nsps)
  var.y <- phylo.var(y, mean.y, phylomat, nsps = nsps)
  
  cor.xy <- as.vector((t(x-mean.x) %*% phylomat%*%(y-mean.x)/(nsps-1))/sqrt(var.x*var.y))
  
  return(cor.xy)
}

calcTraitPairN <- function(data){
  
  vars <- names(data)[!names(data) %in% c("species", "synonyms")]
  
  var.grid <- expand.grid(vars, vars, stringsAsFactors = F)
  var.grid <- var.grid[var.grid[,1] != var.grid[,2],]
  
  indx <- !duplicated(t(apply(var.grid, 1, sort))) # finds non - duplicates in sorted rows
  var.grid <- var.grid[indx, ]
  
  countN <- function(x, data){sum(complete.cases(data[,c(x[1], x[2])]))}
  
  var.grid <- data.frame(var.grid, n = apply(var.grid, 1, FUN = countN, data = data))
  
}

comparePhyloCor <- function(x, data, phylomat, match.dat, tree){
  
  var1 <- unlist(x[1])
  var2 <- unlist(x[2])
  
  data <- data[, c("species", "synonyms", var1, var2)] 
  data <- data[complete.cases(data),]
  
  spp <- data$species
  nsps <- length(spp)
  
  x <- getNamedVector(var1, data = data)
  y <- getNamedVector(var2, data = data)
  
  # Std. correlation
  
  cor <- cor(x, y)
  
  
  #METHOD 1 using phylogenetic relatedness matrix:
  #----------------------------------------------------
  
  #sub.phymat <- subsetPhylomat(spp, phylomat, match.dat)
  
  sub.tree <- drop.tip(tree, tip = tree$tip.label[!tree$tip.label %in% match.dat$synonyms[match(spp, match.dat$species)]])
  sub.phymat <-solve(vcv.phylo(sub.tree))
  spp.m <- match.dat$species[match(sub.tree$tip.label, match.dat$synonyms)]
  x <- x[match(spp.m, names(x))]
  y <- y[match(spp.m, names(y))]
  
  dimnames(sub.phymat) <- list(spp.m, spp.m)
  
  phylocor1 <- getPhyloCor(x, y, phylomat = sub.phymat, nsps = nsps)
  
  
  
  ################################################################################
  
  #METHOD 2 extracting from a PGLS
  #----------------------------------------------------
  
  cd <- comparative.data(phy = tree, data = data, names.col = "synonyms")
  result.pgls <- try(pgls(as.formula(paste(var1, "~", var2, sep = "")), data = cd), silent = T)
  if(class(result.pgls) == "try-error"){phylocor2 <- NA}else{
    
    t <- summary(result.pgls)$coefficients[var2,3]
    df <- as.vector(summary(result.pgls)$fstatistic["dendf"])
    phylocor2 <- sqrt((t*t)/((t*t)+df))*sign(summary(result.pgls)$coefficients[var2,1])
  }
  
  
  return(data.frame(var1 = var1, var2 = var2, cor = cor, phylocor1 = phylocor1, 
                    phylocor2 = phylocor2, n = nsps))
}

pglsPhyloCor <- function(x, data, match.dat, tree, log.vars){
  
  var1 <- unlist(x[1])
  var2 <- unlist(x[2])
  
  data <- data[, c("species", "synonyms", var1, var2)] 
  data <- data[complete.cases(data),]
  
  spp <- data$species
  nsps <- length(spp)
  
  if(var1 %in% log.vars & all(data[,var1] > 0)){
    data[,var1] <- log(data[,var1])
    names(data)[names(data) == var1] <- paste("log", var1, sep = "_")
    var1 <- paste("log", var1, sep = "_")
    
  }
  if(var2 %in% log.vars & all(data[,var2] > 0)){
    data[,var2] <- log(data[,var2])
    names(data)[names(data) == var2] <- paste("log", var2, sep = "_")
    var2 <- paste("log", var2, sep = "_")
  }
  
  # Std. correlation
  cor <- cor(data[,var1], data[,var2])
  
  
  #METHOD 2 extracting from a PGLS
  #----------------------------------------------------
  
  cd <- comparative.data(phy = tree, data = data, names.col = "synonyms", vcv=F)
  
  result.pgls <- try(pgls(as.formula(paste(var1, "~", var2, sep = "")), data = cd, lambda="ML"), 
                     silent = T)
  
  
  if(class(result.pgls) == "try-error"){
    phylocor2 <- NA
    lambda <- NA
    aicc <- NA
  }else{
    
    t <- summary(result.pgls)$coefficients[var2,3]
    df <- as.vector(summary(result.pgls)$fstatistic["dendf"])
    phylocor2 <- sqrt((t*t)/((t*t)+df))*sign(summary(result.pgls)$coefficients[var2,1])
    lambda <- result.pgls$param["lambda"]
    aicc <- result.pgls$aicc
  }
  
  
  return(data.frame(var1 = var1, var2 = var2, cor = cor, phylocor = phylocor2, n = nsps,
                    lambda = lambda, aicc = aicc))
}



