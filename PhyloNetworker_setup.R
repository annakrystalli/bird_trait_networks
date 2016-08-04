# ---- setup_edge ----

# PACKAGES & FUNCTIONS ###############################################################

# packages might need installing
library(shiny)


# SETTINGS ###############################################################

an.ID <- "all"
#an.ID <- "all"

min.n <- 10
log <- T
cutoff <- 0.35

if(log){log.vars <- metadata$code[as.logical(metadata$log)]
}else{log.vars <- ""}


# FILES ##################################################################

wide <- read.csv(file = paste(input.folder,"csv/master wide.csv", sep = ""), fileEncoding = "mac")
spp100 <- unlist(read.csv(file = paste(input.folder,"csv/100spp.csv", sep = "")))

spp.list <- data.frame(species = unique(wide$species))

# trees <- read.tree(file = "tree/Stage2_MayrAll_Hackett_set10_decisive.tre")
# tree <- trees[[1]]
# save(tree, file = "tree/tree.RData")

# LOAD TREE
load(file = paste(input.folder, "tree/tree.RData", sep = ""))

#load match data
load(file = paste(input.folder,"r data/match data/tree m.RData", sep = ""))
phylo.match <- m$data

# WORKFLOW ###############################################################

## DATA ##

### SUBSET TO 100 SPECIES >>>

if(an.ID == "100spp"){
  data <- wide[wide$species %in% spp100,]}else{data <- wide}

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
cat.dat <- m1DataPrep(data = data, datType = c("Cat", "Nom"), 
                      phylo.match = phylo.match)

cat.vg <- varGridGen(cat.dat)





vg <- rbind(data.frame(num.vg, varTypes = "nn"),
            data.frame(bin.vg, varTypes = "bb"),
            data.frame(cat.vg, varTypes = "cc"))


