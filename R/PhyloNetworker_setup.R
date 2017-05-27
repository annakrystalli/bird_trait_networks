# ---- pns-load-files ----
data <- read.csv(file = paste(input.folder,"analytical/master wide.csv", sep = ""))
spp100 <- unlist(read.csv(file = paste(input.folder,"taxo/100spp.csv", sep = "")))
# load selected tree (see select_tree.R)
load(file = paste(input.folder, "tree/tree.RData", sep = ""))
# load match data
load(file = paste(input.folder,"taxo/match data/tree m.RData", sep = ""))
phylo.match <- m$data
# load proc.pars
proc.pars <- read.csv(paste0(input.folder, "metadata/vselect.csv" ), 
                      na.strings = "")
alco.vars <- na.omit(proc.pars$code[proc.pars$allom.cor == T])


# ---- pns-remove-dup-dftips ----
if(remove_dtips){
  phylo.match <- phylo.match[match(data$species, phylo.match$species),]
  dup.syns <- duplicated(phylo.match$synonyms)
  phylo.match <- phylo.match[!dup.syns,]
  data <- data[!dup.syns,]
  tree <- drop.tip(tree, setdiff(tree$tip.label, phylo.match$synonyms))
}

# ---- pns-sex-data-merge ----
if(!exists("sex_dm")){
  sex_dm <- F}
if(sex_dm){
  data <- sex_data_merge(data, proc.pars)
}


# ---- pns-var.delete ----
if(!exists("var.delete")){
  var.delete <- F}
if(var.delete){
  data <- data[, !names(data) %in% proc.pars$code[!proc.pars$phylocor]]
}

# ---- pns-rename-tip.labels ----
tree$tip.label <- phylo.match$species[match(tree$tip.label, phylo.match$synonyms)]

# ---- pns-create-spp.list ----
spp.list <- data.frame(species = unique(data$species))

# ---- pns-subset-data ----
if(an.ID == "100spp"){
  data <- data[data$species %in% spp100,]}
data <- data[!names(data) %in% c("family", "genus", "order")]

# ---- pns-create-analysis-metadata ----
an_meta <- get_an_meta(metadata, data, proc.pars)

# ---- pns-log-data ----
# create log.vars
if(log){log.vars <- an_meta$code[as.logical(an_meta$log)]
# correct-zeros 
valid.log <- validate_log.vars(log.vars, data = data)
zero.vars <- names(valid.log)[!valid.log]
data <- corr_zero.logs(data, replace = zero.vars[!zero.vars %in% "disperse.dist"]) %>%
  corr_zero.logs(replace = "disperse.dist", replace.with = 0.001)
# all log.vars valid to log?
all(validate_log.vars(log.vars, data = data))
# log vars
for(log.var in log.vars){
  data[, log.var] <- log(data[,log.var])
  names(data)[names(data) == log.var] <- paste0(log.var, "_log")
  an_meta$code[an_meta$code == log.var] <- paste0(log.var, "_log")
  if(log.var %in% alco.vars){
    alco.vars[alco.vars == log.var] <- paste0(log.var, "_log")}
}
}else{log.vars <- ""}

# ---- pns-create-ref.obj ----
vg <- varGridGen(data[,!names(data) %in% "species"])
ms_vars <- vg[,1:2] %>% as.matrix() %>% t() %>% as.vector() %>% unique() %>% sort()

# ---- pns-trim-data ----
data <- data[, c("species", ms_vars)]
an_meta <- an_meta %>% filter(code %in% ms_vars)

# ---- pns-correct-cat ----
for(var in c("bite", "eumelan.color.sc", "juv.plum.dist")){
  data[, var] <- round(data[, var])
}
rm("var")

# ---- pns-types-ref ----
vtypes <- data.frame(mgm = c("c", "c", "g", "p", "c"),
                     meta = names(table(
                       an_meta$type[match(ms_vars, an_meta$code)])))
meta_types <- an_meta$type[match(ms_vars, an_meta$code)] %>%
  setNames(ms_vars)
mgm_types <- vtypes$mgm[match(meta_types, vtypes$meta)] %>%
  setNames(ms_vars)


# ---- pns-center_dt ----
center_vars <- which(names(data) %in% names(mgm_types)[mgm_types != "c"])
data[,center_vars] <-  scale(data[,center_vars], center = T, scale = F)

# ---- pns-create-ref.obj2 ----
ms_spp <- data$species
row.names(data) <- ms_spp
spp.ranks <- getSppRanks(ms_vars, data = data, load = T, input.folder, an.ID)
