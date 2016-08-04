# SETUP ###############################################################

rm(list=ls())

# DOWNLOAD FILE FROM GITHUB, UPDATE AND SOURCE FROM APPROPRIATE PATH
source("~/Documents/workflows/bird_trait_networks/setup.R")
setwd(input.folder) # e.g. googledrive/bird trait networks/inputs/data. Set in setup script

# PACKAGES & FUNCTIONS ###############################################################

source("~/Documents/workflows/bird_trait_networks/PhyloNetworker_functions.R")

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






# PHYLOGENTICALLY CORRECTED CORRELATIONS ####################################################
num.res <- NULL
for(i in 1:dim(num.vg)[1]){
  
  num.res <- rbind(num.res, pglsPhyloCor(x = num.vg[i, 1:2], data = num.dat, 
                                         tree = tree, log.vars = log.vars, 
                                         datTypes = "nn"))
  print(i)
}


num.res <- num.res[order(abs(num.res$phylocor), decreasing = T),]

write.csv(num.res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_num.csv", sep = ""),
          row.names = F)


bin.res <- NULL
for(i in 1:dim(bin.vg)[1]){
  
  bin.res <- rbind(bin.res, pglsPhyloCor(x = bin.vg[i, 1:2], data = bin.dat, 
                                         tree = tree, log.vars = log.vars, 
                                         datTypes = "bb"))
  print(i)
}


bin.res <- bin.res[order(abs(bin.res$phylocor), decreasing = T),]

write.csv(bin.res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_bin.csv", sep = ""),
          row.names = F)

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

