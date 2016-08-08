# ---- get-log.vars ----
if(log){log.vars <- metadata$code[as.logical(metadata$log)]
}else{log.vars <- ""}

# ---- load-files
wide <- read.csv(file = paste(input.folder,"csv/master wide.csv", sep = ""), fileEncoding = "mac")
spp100 <- unlist(read.csv(file = paste(input.folder,"csv/100spp.csv", sep = "")))
spp.list <- data.frame(species = unique(wide$species))
# load selected tree (see select_tree.R)
load(file = paste(input.folder, "tree/tree.RData", sep = ""))
# load match data
load(file = paste(input.folder,"r data/match data/tree m.RData", sep = ""))
phylo.match <- m$data

# ---- subset-data ----
if(an.ID == "100spp"){
  data <- wide[wide$species %in% spp100,]}else{data <- wide}

# ---- num-prep ---- 
num.dat <- m1DataPrep(data = data, datType = c("Int", "Con"), 
                      phylo.match = phylo.match)


num.vg <- varGridGen(num.dat)
## make sure variables to be logged are > 0
log.vars <- log.vars[sapply(log.vars, FUN = function(x, dat){all(na.omit(dat[,x]) > 0)},
                            dat = num.dat)]
# ---- bin-prep ---- 
bin.dat <- m1DataPrep(data = data, datType = "Bin", 
                      phylo.match = phylo.match)
bin.vg <- varGridGen(bin.dat)


# ---- cat-prep ----
cat.dat <- m1DataPrep(data = data, datType = c("Cat", "Nom"), 
                      phylo.match = phylo.match)
cat.vg <- varGridGen(cat.dat)




# ---- full.vg ----
vg <- rbind(data.frame(num.vg, varTypes = "nn"),
            data.frame(bin.vg, varTypes = "bb"),
            data.frame(cat.vg, varTypes = "cc"))


