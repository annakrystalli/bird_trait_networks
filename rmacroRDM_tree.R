## ----global-setup, echo = F----------------------------------------------
rm(list = ls(all.names = T))
options(stringsAsFactors = F)
.rs.restartR()

if(exists("file_setup_path")){}else{
  file_setup_path <- "file_setup.R"
  source(file_setup_path)}

## ----intialise-project ----
wkf = "rmacro"
source(paste0(script.folder, "project_ui.R"))

## ----master-configuration, eval=T----------------------------------------
init_db(data.folder = "/Users/Anna/Google Drive/bird trait networks/",
        script.folder = "~/Documents/workflows/bird_trait_networks/", 
        spp.list_src = "D0")

# SETTINGS ###############################################################


# match phylogenetic tree ##################################################################

# ---- pm-spp.list-from-master ----
load(file = paste(output.folder, "data/master.rda", sep = ""))
spp.list <- master$spp.list
syn_links <- attr(spp.list, "syn_links")

load_sys.ref(view = F)



# trees <- read.tree(file = "tree/Stage2_MayrAll_Hackett_set10_decisive.tre")
# tree <- trees[[1]]
# save(tree, file = "tree/tree.RData")
load(file = paste0(input.folder, "tree/tree.RData"))

# WORKFLOW ###############################################################

## Match tree species names to master species list
treespp <- data.frame(species = tree$tip.label)


# CREATE MASTER


m <- matchObj(dcode = "tree", spp.list = spp.list,
              data = treespp,
              sub = "spp.list", 
              meta = NULL) # use addMeta function to manually add metadata.

m$meta$ref <- "hackett"

m <- m %>% 
  separateDatMeta() %>% 
  compileMeta() %>%
  checkVarMeta() %>%
  dataMatchPrep()

## ----data-spp-match------------------------------------------------------
m <- dataSppMatch(m, addSpp = T)
save(m, file = paste(output.folder, "data/match objects/", m$dcode, "m.rda", sep = ""))


m <- testSynonym("Buphagus erythrorhynchus", m)
syn_links <- updateSynlinks(syn_links, m)

write.csv(syn_links, "taxo/syn_links.csv", row.names = F)


save(m, file = "r data/match data/tree m.RData")

################################################################################

#load match data
load(file = "r data/match data/tree m.RData")
match.dat <- m$data
