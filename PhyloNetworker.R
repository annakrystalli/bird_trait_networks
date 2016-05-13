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

calcTraitPairN <- function(data){
  
  vars <- names(data)[!names(data) %in% c("species", "synonyms")]
  
  var.grid <- expand.grid(vars, vars, stringsAsFactors = F)
  var.grid <- var.grid[var.grid[,1] != var.grid[,2],]
  
  indx <- !duplicated(t(apply(var.grid, 1, sort))) # finds non - duplicates in sorted rows
  var.grid <- var.grid[indx, ]
  
  countN <- function(x, data){sum(complete.cases(data[,c(x[1], x[2])]))}
  
  var.grid <- data.frame(var.grid, n = apply(var.grid, 1, FUN = countN, data = data))
  
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

# SETTINGS ###############################################################
dir.create(paste(output.folder, "data/phylocors/", sep = ""))
dir.create(paste(output.folder, "data/networks/", sep = ""))

an.ID <- "100spp"
min.n <- 10
log <- T
if(log){log.vars <- metadata$code[as.logical(metadata$log)]}else{log.vars <- ""}


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
match.dat <- m$data

# WORKFLOW ###############################################################

## NUMERIC VARIABLES >>>>


# separate numeric variables
num.var <- metadata$code[metadata$type %in% c("Int", "Con")]
num.dat <- wide[,c("species", names(wide)[names(wide) %in% num.var])]

#Remove duplicate species matching to the same species on the tree
num.dat <- num.dat[num.dat$species %in% match.dat$species[match.dat$data.status != "duplicate"],]

# add synonym column to data 
num.dat$synonyms <- match.dat$species[match(num.dat$species, match.dat$species)]


## SUBSET TO 100 SPECIES >>>

if(an.ID == "100spp"){
  num.dat <- num.dat[num.dat$species %in% spp100,]}


# VARIABLES COMBINATION DATA AVAILABILITY >>>

## Create grid of unique variable combinations, calculate data availability for each and sort


var.grid <- calcTraitPairN(num.dat)
var.grid <- var.grid[var.grid$n > min.n,]
var.grid <- var.grid[order(var.grid$n, decreasing = T),]


# PHYLOGENTICALLY CORRECTED CORRELATIONS ####################################################


res <- NULL
for(i in 1:dim(var.grid)[1]){
  
  res <- rbind(res, pglsPhyloCor(var.grid[i, 1:2], data = num.dat, 
                                 match.dat = match.dat, tree = tree, log.vars = log.vars))
  print(i)
}


res <- res[order(abs(res$phylocor), decreasing = T),]

write.csv(res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                     min.n, if(log){"_log"},".csv", sep = ""),
          row.names = F)



## NETWORK ANALYSIS ####################################################################

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




# igraph ###############################################################################

load(file = paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                       if(log){"_log"},".RData", sep = ""))


edges <- res[abs(res$phylocor) > 0.35 & !is.na(res$phylocor),]

g.n <- as.vector(apply(edges[,c("var1", "var2")],1,FUN = unlist))
n.names <- unique(g.n)
g.code <- match(g.n, n.names)  

# node colours
role.levels <- levels(unique(net[[1]]$role))
role.codes <- 
role.cols <- RColorBrewer::brewer.pal(7, "BrBG")
n.cols <- role.cols[match(as.character(net[[1]]$role)[match(n.names, net[[1]]$name)], role.levels)]

# edge colours
edge.cols <- rep("black", length(edges$phylocor))
edge.cols[edges$phylocor >= 0.1] <- "red"
edge.cols[edges$phylocor < -0.1] <- "blue"


require(igraph)

G <- graph(g.code, directed = FALSE )

# Assign attributes to the graph
G$name    <- "correlation network of numeric bird traits"

# Assign attributes to the graph's vertices
V(G)$name  <- 1:length(n.cols)
V(G)$color <- n.cols
E(G)$edge.color <- edge.cols


# Assign attributes to the edges
E(G)$weight <- edges$phylocor



# Plot the graph -- details in the "Drawing graphs" section of the igraph manual


png(filename = paste(output.folder, "figures/net.png", sep = ""),
    width =  1500, height = 1194, pointsize = 20)

par(mai = c(0,0.5,0.5,0))

plot(G, layout = layout.fruchterman.reingold, 
     main = G$name,
     vertex.label = V(G)$name,
     vertex.size = 5.5,
     vertex.color= V(G)$color,
     vertex.frame.color= "white",
     vertex.label.color = "black",
     vertex.label.family = "sans",
     edge.width=E(G)$weight,
     vertex.label.cex=0.65,
     edge.color= edge.cols)



legend("bottomleft", legend = role.levels[which(role.levels %in% as.character(net[[1]]$role))], 
       fill = role.cols[which(role.levels %in% as.character(net[[1]]$role))], 
       bty = "n", cex = 1)

dev.off()


       
#source('http://bioconductor.org/biocLite.R')
#biocLite ('RCytoscape')

