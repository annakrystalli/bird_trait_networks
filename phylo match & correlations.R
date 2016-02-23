# SETUP ###############################################################
rms <- ls()[ls() != "trees"]
rm(list=rms)
source("~/Documents/workflows/bird_trait_networks/setup.R")
setwd(input.folder)

# PACKAGES & FUNCTIONS ###############################################################

library(caper)
library(geiger)
library(PHYLOGR)

source(paste(script.folder, "functions.R", sep = ""))

# SETTINGS ###############################################################

qcmnames = c("qc", "observer", "ref", "n", "notes")
taxo.var <- c("species", "order","family", "subspp", "parent.spp")
var.var <- c("var", "value", "data")
var.omit <- c("no_sex_maturity_d", "adult_svl_cm", "male_maturity_d")


# FILES ##################################################################

#trees <- read.tree(file = "tree/Stage2_MayrAll_Hackett_set10_decisive.tre")
D0 <- read.csv(file = "csv/D0.csv" ,fileEncoding = "mac")


D0 <- D0[!D0$species == "Conuropsis_carolinensis",] #species extinct
D0 <- D0[!D0$species == "Podilymbus_gigas",] #species extinct
D0 <- D0[!D0$species == "Xenicus_longipes",] #species extinct
D0 <- D0[!D0$species == "Ectopistes_migratorius",] #species extinct
D0 <- D0[!D0$species == "Pezophaps_solitaria",] #species extinct
D0 <- D0[!D0$species == "Pinguinus_impennis",] #species extinct
D0 <- D0[!D0$species == "Porphyrio_albus",] #species extinct


spp.list <- data.frame(species = unique(D0$species))
synonyms  <- read.csv("r data/synonyms.csv", stringsAsFactors=FALSE)
r.synonyms <- read.csv("r data/bird_species_names.csv", stringsAsFactors=FALSE)
m.synonyms <- read.csv("r data/match data/tree mmatched.csv", stringsAsFactors=FALSE)


# WORKFLOW ###############################################################

## Match tree species names

tree <- trees[[1]]
treespp <- data.frame(species = tree$tip.label)

dl <- processDat(file = NULL, dat = treespp, label = F, taxo.dat, var.omit, input.folder,
           observer = NULL, qc = NULL, ref = "Hackett", n = NULL, notes = NULL,
           master.vname = "master.vname")

m <- matchObj(data.ID = "tree", spp.list = spp.list, data = dl$data, 
              status = "unmatched", 
                          sub = "spp.list",
                          qcref = dl$qcref)

unmatched <- m$spp.list$species[!m$spp.list$species %in% m$data$species]

m <- dataSppMatch(m, unmatched, ignore.unmatched = T, synonyms = r.synonyms, trim.dat = F, 
                  retain.dup = F)

unmatched <- m$spp.list$species[!m$spp.list$species %in% m$data$species]

m <- dataSppMatch(m, unmatched, ignore.unmatched = T, synonyms = synonyms, 
                  trim.dat = F, retain.dup = F)


unmatched <- m$spp.list$species[!m$spp.list$species %in% m$data$species]

m <- dataSppMatch(m, unmatched, ignore.unmatched = F, synonyms = m.synonyms, 
                  trim.dat = F, retain.dup = F)

################################################################################

#Prepare phylogenetic relatedness matrix

tdat <- m$data


invC <-solve(vcv.phylo(tree))


## vars 

var.n <- sort(table(D0$var), decreasing = T)
var1 <- names(var.n)[1]
var2 <- names(var.n)[2]

x.dat <- D0[D0$var == var1,]
y.dat <- D0[D0$var == var2,]

x <- as.numeric(x.dat$value[match(spp, x.dat$species)])
y <- as.numeric(x.dat$value[match(spp, y.dat$species)])


spp <- intersect(unique(D0$species[D0$var == var1]), unique(D0$species[D0$var == var2]))
nsps <- length(spp)


tdat$synonyms[match(spp, tdat$species)]

vmat <- matrix(NA, nrow = nsps, ncol = nsps, dimnames = list(spp,spp))


m.id <- cbind(rep(spp, times = nsps), rep(spp, each = nsps))
p.id <- cbind(rep(tdat$synonyms[match(spp, tdat$species)], times = nsps), 
              rep(tdat$synonyms[match(spp, tdat$species)], each = nsps))

# Produce new phylo.matrix for species available for variable combination
vmat[m.id] <- invC[p.id]



mean.x <-colSums(vmat%*%x)/sum(vmat)
mean.y <-colSums(vmat%*%y)/sum(vmat)

var.x <-t(x-mean.x) %*% vmat%*%(x-mean.x)/(nsps-1)
var.y <-t(y-mean.x) %*% vmat%*%(y-mean.x)/(nsps-1)
cor.xy <-(t(x-mean.x) %*% vmat%*%(y-mean.x)/(nsps-1))/sqrt(var.x*var.y)
cor.xy




#save(tree, file = paste(input.folder, "tree/example tree.RData", sep = ""))

# rTraitCont: simulates the evolution of a continuous character along a phylogeny
# rescale: applies transformations to a phylo tree. lambda = is one of the Pagel (1999) 
# models that fits the extent to which the phylogeny predicts covariance among trait values for species. The model effectively transforms the tree as follows: values of lambda near 0 cause the phylogeny to become more star-like, and a lambda value of 1 recovers the BM model. The parameter used for transformation is lambda.

x.pic<- pic(x, tree)
y.pic<- pic(y, tree)
cor.origin(x.pic,y.pic) # this will equal the cor.xy below


mean.x <-colSums(invC%*%x)/sum(invC)
mean.y <-colSums(invC%*%y)/sum(invC)
vector.ones<-as.matrix(rep(1,nsps))
var.x <-t(x-vector.ones%*%mean.x) %*% invC%*%(x-vector.ones%*%mean.x)/(nsps-1)
var.y <-t(y-vector.ones%*%mean.x) %*% invC%*%(y-vector.ones%*%mean.x)/(nsps-1)
cor.xy <-(t(x-vector.ones%*%mean.x) %*% invC%*%(y-vector.ones%*%mean.x)/(nsps-1))/sqrt(var.x*var.y)
cor.xy

