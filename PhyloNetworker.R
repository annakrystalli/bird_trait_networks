# SETUP ###############################################################

rm(list=ls())

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

calcTraitPairN <- function(data){
  
  vars <- names(data)[!names(data) %in% c("species", "synonyms")]
  
  var.grid <- expand.grid(vars, vars, stringsAsFactors = F)
  var.grid <- var.grid[var.grid[,1] != var.grid[,2],]
  
  indx <- !duplicated(t(apply(var.grid, 1, sort))) # finds non - duplicates in sorted rows
  var.grid <- var.grid[indx, ]
  
  countN <- function(x, data){sum(complete.cases(data[,c(x[1], x[2])]))}
  
  var.grid <- data.frame(var.grid, n = apply(var.grid, 1, FUN = countN, data = data))
  
}

pglsPhyloCor <- function(x, data, phylo.match, tree, log.vars){
  
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
  
  result.pgls <- try(pgls(as.formula(paste(var1, "~", var2, sep = "")), data = cd, lambda="ML"))
  
  
  if(class(result.pgls) == "try-error"){
    
    phylocor2 <- NA
    lambda <- NA
    error <- gsub(pattern = ".*\n  ", "", geterrmessage())
    
  }else{
    
    t <- summary(result.pgls)$coefficients[var2,3]
    df <- as.vector(summary(result.pgls)$fstatistic["dendf"])
    phylocor2 <- sqrt((t*t)/((t*t)+df))*sign(summary(result.pgls)$coefficients[var2,1])
    lambda <- result.pgls$param["lambda"]
    error <- NA
    
  }
  
  
  return(data.frame(var1 = var1, var2 = var2, cor = cor, phylocor = phylocor2, n = nsps,
                    lambda = lambda, error = error))
}

m1DataPrep <- function(data, datType, phylo.match, ...) {
  
  var <- metadata$code[metadata$type %in% datType]
  dat <- data[,c("species", names(data)[names(data) %in% var])]
  
  #Remove duplicate species matching to the same species on the tree
  dat <- dat[dat$species %in% phylo.match$species[phylo.match$data.status != "duplicate"],]
  
  # add synonym column to data 
  dat$synonyms <- phylo.match$species[match(dat$species, phylo.match$species)]
  
  return(dat)
  
}

## Create grid of unique variable combinations, calculate data availability for
## each and sort. min.n needs to be set in working environment.
varGridGen <- function(dat, ...) {
  
  var.grid <- calcTraitPairN(dat)
  var.grid <- var.grid[var.grid$n > min.n,]
  var.grid <- var.grid[order(var.grid$n, decreasing = T),] 
  
  return(var.grid)
  
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


## Extracting interaction weights

### Extracting correlations >>>>

#### numeric variables 

num.dat <- m1DataPrep(data = wide, datType = c("Int", "Con"), 
                      phylo.match = phylo.match)

# separate numeric variables



## SUBSET TO 100 SPECIES >>>

if(an.ID == "100spp"){
  num.dat <- num.dat[num.dat$species %in% spp100,]}


# VARIABLES COMBINATION DATA AVAILABILITY >>>

num.vg <- varGridGen(num.dat)

## make sure variables to be logged are > 0
log.vars <- log.vars[sapply(log.vars, FUN = function(x, dat){all(na.omit(dat[,x]) > 0)},
                             dat = num.dat)]

# PHYLOGENTICALLY CORRECTED CORRELATIONS ####################################################


res <- NULL
for(i in 1:dim(var.grid)[1]){
  
  res <- rbind(res, pglsPhyloCor(var.grid[i, 1:2], data = num.dat, 
                                 phylo.match = phylo.match, tree = tree, log.vars = log.vars))
  print(i)
}


res <- res[order(abs(res$phylocor), decreasing = T),]

write.csv(res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                     min.n, if(log){"_log"},".csv", sep = ""),
          row.names = F)


#### binary variables

# separate numeric variables
bin.var <- metadata$code[metadata$type == "Bin"]
bin.dat <- wide[,c("species", names(wide)[names(wide) %in% bin.var])]

var.grid <- calcTraitPairN(bin.dat)
var.grid <- var.grid[var.grid$n > min.n,]
var.grid <- var.grid[order(var.grid$n, decreasing = T),]

fitPagel(tree, x, y, method="fitMk", ...)


### mixed models





#################################################################################
## NETWORK ANALYSIS ####################################################################
##################################################################################
res <- res[abs(res$phylocor) > cutoff & !is.na(res$phylocor),]

library(rnetcarto)
net.d <- res[!is.na(res$phylocor), c("var1", "var2", "phylocor")]
net.list <- as.list(net.d)
net <- netcarto(web = net.list, seed = 1)


write.csv(cbind(net[[1]], modularity = net[[2]]), 
          paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                if(log){"_log"},".csv", sep = ""),
          row.names = F)

save(net, file = paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                       if(log){"_log"},".RData", sep = ""))






       
#source('http://bioconductor.org/biocLite.R')
#biocLite ('RCytoscape')

