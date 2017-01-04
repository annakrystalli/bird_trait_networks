# ---- pns-get-log.vars ----
if(log){log.vars <- metadata$code[as.logical(metadata$log)]
}else{log.vars <- ""}

# ---- pns-load-files ----
wide <- read.csv(file = paste(input.folder,"csv/master wide.csv", sep = ""), fileEncoding = "mac")
spp100 <- unlist(read.csv(file = paste(input.folder,"csv/100spp.csv", sep = "")))
# load selected tree (see select_tree.R)
load(file = paste(input.folder, "tree/tree.RData", sep = ""))
# load match data
load(file = paste(input.folder,"r data/match data/tree m.RData", sep = ""))
phylo.match <- m$data


# ---- pns-remove-dup-dftips ----
if(remove_dtips){
  phylo.match <- phylo.match[match(wide$species, phylo.match$species),]
  dup.syns <- duplicated(phylo.match$synonyms)
  phylo.match <- phylo.match[!dup.syns,]
  wide <- wide[!dup.syns,]
  tree <- drop.tip(tree, setdiff(tree$tip.label, phylo.match$synonyms))
}

# ---- pns-rename-tip.labels ----
tree$tip.label <- phylo.match$species[match(tree$tip.label, phylo.match$synonyms)]

# ---- pns-create-spp.list ----
spp.list <- data.frame(species = unique(wide$species))

# ---- pns-subset-data ----
if(an.ID == "100spp"){
  data <- wide[wide$species %in% spp100,]}else{data <- wide}
data <- data[!names(data) %in% c("family", "genus", "order")]

# ---- pns-create-ref.obj ----
vg <- varGridGen(data[,!names(data) %in% "species"])
ms_vars <- vg[,1:2] %>% as.matrix() %>% t() %>% as.vector() %>% unique()

# ---- pns-types-ref ----
vtypes <- data.frame(mgm = c("c", "c", "g", "p", "c"),
                     meta = names(table(
                       metadata$type[match(ms_vars, metadata$code)])))
meta_types <- metadata$type[match(ms_vars, metadata$code)] %>%
  setNames(ms_vars)
mgm_types <- vtypes$mgm[match(meta_types, vtypes$meta)] %>%
  setNames(ms_vars)

# ---- pns-get-TD ----
if(!file.exists(paste(input.folder, "r data/TD_", an.ID, 
                      ".Rdata", sep = ""))){
  td <- getTD(data, tree, all.vars = F)
  TDdf <- extractTD(td, vars = ms_vars)
  save(td, TDdf, file = paste(input.folder, "r data/TD_", an.ID, 
                              ".Rdata", sep = ""))}else{
                                load(file = paste(input.folder, "r data/TD_", 
                                                  an.ID, ".Rdata", sep = ""))}
# ---- pns-create-ref.obj2 ----
ms_spp <- data$species
row.names(data) <- ms_spp
spp.ranks <- getSppRanks(ms_vars, data = data, load = T, input.folder, an.ID)
