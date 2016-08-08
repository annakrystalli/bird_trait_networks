# DOWNLOAD FILE FROM GITHUB, UPDATE AND SOURCE FROM APPROPRIATE PATH
source("~/Documents/workflows/bird_trait_networks/setup.R")
setwd(input.folder) # e.g. googledrive/bird trait networks/inputs/data. Set in setup script

# PACKAGES & FUNCTIONS ###############################################################

# packages might need installing
library(caper)
library(geiger)
library(PHYLOGR)
library(rnetcarto)
library(igraph)
library(phytools)
library(vegan)
library(MCMCglmm)

calcTraitPairN <- function(data){
  
  vars <- names(data)[!names(data) %in% c("species", "synonyms")]
  
  var.grid <- expand.grid(vars, vars, stringsAsFactors = F)
  var.grid <- var.grid[var.grid[,1] != var.grid[,2],]
  
  indx <- !duplicated(t(apply(var.grid, 1, sort))) # finds non - duplicates in sorted rows
  var.grid <- var.grid[indx, ]
  
  countN <- function(x, data){sum(complete.cases(data[,c(x[1], x[2])]))}
  
  var.grid <- data.frame(var.grid, n = apply(var.grid, 1, FUN = countN, data = data))
  
}

pglsPhyloCor <- function(x, data, tree, log.vars = NULL, datTypes){
  
  var1 <- unlist(x[1])
  var2 <- unlist(x[2])
  data <- data[, c("species", "synonyms", var1, var2)]
  TD <- getTD(data, tree)
  data <- data[complete.cases(data),]
  
  
  if(datTypes == "nn") {
    
    if(validateLog(var = "var1", log.vars, data, var1 = var1)){
      data <- logData(data, "var1", var1, var2)
      var1 <- paste("log", var1, sep = "_")
    }
    if(validateLog(var = "var2", log.vars, data, var2 = var2)){
      data <- logData(data, "var2", var1, var2)
      var2 <- paste("log", var2, sep = "_")
    }
    row <- fitNumDat(data, tree, var1, var2, TD)
  }
  
  if(datTypes == "bb") {
    row <-  fitBinDat(data, tree, var1, var2, TD) 
  }
  
  if(datTypes == "cc") {
    
    getLevels <- function(metadata, var) {
      meta <- list(levels = metadata[metadata$code == var, "scores"], 
                   labels = metadata[metadata$code == var, "levels"]) %>%
        sapply(FUN = function(x){strsplit(x, ";")}, simplify = T) %>% 
        data.frame(stringsAsFactors = F)
      
      return(meta)
    }
    
    
    meta <- list(getLevels(metadata, var1), getLevels(metadata, var2))
    
    data[,c(var1, var2)] <- mapply(FUN = function(x, meta){factor(as.character(x),
                                                                  levels = meta$levels,
                                                                  labels = meta$labels)},
                                   x = data[,c(var1, var2)], meta = meta) %>%
      as.data.frame()
    
    row <-  fitCatDat(data, tree, var1, var2, TD) }
  
  return(row)
  
}

fitNumDat <- function(data, tree, var1, var2, TD) {
  
  
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$synonyms))
  tree$tip.label <- data$species[match(tree$tip.label, data$synonyms)] 
  
  physig.var1 <- phylosig(tree, setNames(data[,var1], data$species))
  physig.var2 <- phylosig(tree, setNames(data[,var2], data$species))
  
  
  # Std. correlation
  cor <- cor(data[,var1], data[,var2])
  
  #METHOD extracting from a PGLS
  #----------------------------------------------------
  cd <- comparative.data(phy = tree, data = data, names.col = "synonyms", vcv=F)
  
  result.pgls <- try(pgls(as.formula(paste(var1, "~", var2, sep = "")),
                          data = cd, lambda="ML"))
  
  mod.l <- pgls.profile(result.pgls, 'lambda')
  plot(mod.l)
  plot(result.pgls)
  
  
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

fitBinDat <- function(data, tree, var1, var2, TD) {
  
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$synonyms))
  tree$tip.label <- data$species[match(tree$tip.label, data$synonyms)]  
  
  # Std. correlation
  cor <- NA
  
  data <- data[match(tree$tip.label, data$species),]
  
  x <- setNames(factor(data[,var1]), data$species)
  y <- setNames(factor(data[,var2]), data$species)
  
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
  
  
  physig.var1 <- eval(parse(
    text = paste(
      "phylo.d(data = data, phy = tree, names.col = species, binvar =", var1,")")))$DEstimate
  physig.var2 <- eval(parse(
    text = paste(
      "phylo.d(data = data, phy = tree, names.col = species, binvar =", var2,")")))$DEstimate
  
  
  result.pgl <- try(fitPagel(tree = tree, x = x, y = y, method="fitMk"), silent = T)
  
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

fitCatDat <- function(data, tree, var1, var2, TD) {
  
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
  
  results.mcmc <- MCMCglmm(as.formula(paste(var1, "~ trait + trait:", var2, "-1")), 
                           random = ~us(trait):species,
                           rcov = ~us(trait):units,
                           data = data,
                           ginverse = list(species=Ainv),
                           family = "categorical",
                           prior=prior.phyl,
                           thin   = 50,
                           burnin = 3000,
                           nitt   = 100000)
  
  par(mfrow=c(8,2), mar=c(2,2,1,0))
  plot(results.mcmc$Sol, auto.layout=F)
  
}


getTD <- function(data, tree) {
  
  ptree <- drop.tip(tree, setdiff(tree$tip.label, data$synonyms))
  ptree$tip.label <- data$species[match(ptree$tip.label, data$synonyms)] 
  dist <- cophenetic(ptree)
  dspp <- rownames(dist)
  comm.spp <- data$species[complete.cases(data)]
  comm <- setNames(data.frame(t(dspp %in% comm.spp)), dspp)
  TD <-  taxondive(comm, dist)
  return(TD)
}


m1DataPrep <- function(data, datType, phylo.match, ...) {
  
  var <- metadata$code[metadata$type %in% datType]
  dat <- data[,c("species", names(data)[names(data) %in% var])]
  
  #Remove duplicate species matching to the same species on the tree
  dat <- dat[dat$species %in% phylo.match$species[phylo.match$data.status != "duplicate"],]
  
  # add synonym column to data 
  dat$synonyms <- phylo.match$synonyms[match(dat$species, phylo.match$species)]
  
  return(dat)
  
}

## Create grid of unique variable combinations, calculate data availability for
## each and sort. min.n needs to be set in working environment.
varGridGen <- function(dat, ...) {
  
  var.grid <- calcTraitPairN(dat)
  var.grid <- var.grid[var.grid$n >= min.n,]
  var.grid <- var.grid[order(var.grid$n, decreasing = T),] 
  
  return(var.grid)
  
}

validateLog <- function(var, log.vars, data, var1 = NULL, var2 = NULL) {
  get(var) %in% log.vars & all(data[,get(var)] > 0)
}

logData <- function(data, var, var1, var2) {
  data[,get(var)] <- log(data[,get(var)])
  names(data)[names(data) == get(var)] <- paste("log", get(var), sep = "_")
  return(data)
}

# SETTINGS ###############################################################
dir.create(paste(output.folder, "data/phylocors/", sep = ""))
dir.create(paste(output.folder, "data/networks/", sep = ""))

an.ID <- "100spp"
min.n <- 10
log <- T
cutoff <- 0.35

if(log){log.vars <- metadata$code[as.logical(metadata$log)]
}else{log.vars <- ""}


# FILES ##################################################################

wide <- read.csv(file ="csv/master wide.csv", fileEncoding = "mac")
spp100 <- unlist(read.csv(file ="csv/100spp.csv"))

spp.list <- data.frame(species = unique(wide$species))

# trees <- read.tree(file = "tree/Stage2_MayrAll_Hackett_set10_decisive.tre")
# tree <- trees[[1]]
# save(tree, file = "tree/tree.RData")

# LOAD TREE
load(file = "tree/tree.RData")

#load match data
load(file = "r data/match data/tree m.RData")
phylo.match <- m$data

# WORKFLOW ###############################################################

## DATA ##

### SUBSET TO 100 SPECIES >>>

if(an.ID == "100spp"){
  data <- wide[wide$species %in% spp100,]}

#### numeric variables 

num.dat <- m1DataPrep(data = data, datType = c("Int", "Con"), 
                      phylo.match = phylo.match)
num.vg <- varGridGen(num.dat)

## make sure variables to be logged are > 0
log.vars <- log.vars[sapply(log.vars, FUN = function(x, dat){all(na.omit(dat[,x]) > 0)},
                            dat = num.dat)]


#### binary variables
bin.dat <- m1DataPrep(data = data, datType = "Bin", 
                      phylo.match = phylo.match)

bin.vg <- varGridGen(bin.dat)


#### categorical variables
cat.dat <- m1DataPrep(data = data, datType = "Cat", 
                      phylo.match = phylo.match)

cat.vg <- varGridGen(cat.dat)





vg <- rbind(data.frame(num.vg, varTypes = "nn"),
            data.frame(bin.vg, varTypes = "bb"),
            data.frame(cat.vg, varTypes = "cc"))


