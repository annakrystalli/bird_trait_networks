## ----global-setup, echo = F----------------------------------------------
rm(list=ls())
options(stringsAsFactors = F)

if(exists("file_setup_path")){}else{
  file_setup_path <- "file_setup.R"
  source(file_setup_path)}
source(paste0(script.folder,"project_ui.R"))

## ----intialise-project ----
wkf = "rmacro"
param = "rmacro.R"
source("~/Documents/workflows/bird_trait_networks/project_ui.R")

## ----master-configuration, eval=T----------------------------------------
init_db(data.folder = "/Users/Anna/Google Drive/bird trait networks/",
        script.folder = "~/Documents/workflows/bird_trait_networks/", 
        spp.list_src = "D0", fileEncoding = "mac")

# SETTINGS ###############################################################


# match phylogenetic tree ##################################################################

# ---- pm-spp.list-from-master ----
load(file = paste(output.folder, "data/master.RData", sep = ""))
spp.list <- master$spp.list

fcodes <- ensure_fcodes(meta.vars)
load_sys.ref(view = T)
syn_links <- read.csv(paste0(input.folder,"taxo/syn_links.csv"), stringsAsFactors = F)


# trees <- read.tree(file = "tree/Stage2_MayrAll_Hackett_set10_decisive.tre")
# tree <- trees[[1]]
# save(tree, file = "tree/tree.RData")
load(file = "tree/tree.RData")

# WORKFLOW ###############################################################

## Match tree species names to master species list
treespp <- data.frame(species = tree$tip.label)


# CREATE MASTER


m <- matchObj(data.ID = "tree", spp.list = spp.list,
              data = treespp,
              sub = "spp.list", filename = NULL, 
              meta = NULL) # use addMeta function to manually add metadata.

m <- processDat(m, input.folder, var.omit)
m$meta$ref <- "hackett"
m <-  checkVarMeta(m, master$metadata) %>%  dataMatchPrep()
m <- dataSppMatch(m, syn_links = syn_links, addSpp = F, ignore.unmatched = F)


m <- testSynonym("Buphagus erythrorhynchus", m)
syn_links <- updateSynlinks(syn_links, m)

write.csv(syn_links, "taxo/syn_links.csv", row.names = F)




# Match data set to spp.list and process
output <- masterDataFormat(m, meta.vars, match.vars, var.vars)

write.csv(output$data, file =  "csv/D1.long.csv", row.names = F, fileEncoding = "mac")

spp.list <- output$spp.list


save(m, file = "r data/match data/tree m.RData")

################################################################################

#load match data
load(file = "r data/match data/tree m.RData")
match.dat <- m$data

# Create wide dataset:
wide <- widenMaster(vars = unique(D0$var), species = unique(D0$species), 
                    master = D0, metadata = metadata)

# separate numeric variables
num.var <- metadata$master.vname[metadata$type %in% c("Int", "Con")]
num.dat <- wide[,c("species", names(wide)[names(wide) %in% num.var])]

#Remove duplicate species matching to the same species on the tree
num.dat <- num.dat[num.dat$species %in% match.dat$species[match.dat$data.status != "duplicate"],]

# add synonym column to data 
num.dat$synonyms <- match.dat$species[match(num.dat$species, match.dat$species)]

# VARIABLES

## Create grid of unique variable combinations, calculate data availability for each and sort
var.grid <- calcTraitPairN(num.dat)
var.grid <- var.grid[var.grid$n > 30,]
var.grid <- var.grid[order(var.grid$n, decreasing = T),]

#Prepare phylogenetic relatedness matrix

# phylomat <-solve(vcv.phylo(tree))  #REALLY TIME CONSUMING
# save(phylomat, file = "tree/phylomat.RData")
# load(file = "tree/phylomat.RData")

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

res <- NULL
for(i in 1:dim(res.grid)[1]){
  
  res <- rbind(res, comparePhyloCor(res.grid[i,1:2], data = num.dat, 
               phylomat = phylomat, match.dat = match.dat, tree = tree))
}


write.csv(res, "r data/phylocor comps.csv", row.names = F)


res2 <- apply(var.grid[1:5,], 2, FUN = comparePhyloCor, data = num.dat, 
                phylomat = phylomat, match.dat = match.dat)



res.i <- read.csv(paste(input.folder, "r data/phylocor comps.csv", sep = ""), stringsAsFactors = F)
res.grid <- res.i[order(abs(res.i$phylocor1 - res.i$phylocor2), decreasing = T),][1:400,c(1,2,6)]


dim(var.grid)[1]









