## ----global-setup, echo = F----------------------------------------------
rm(list=ls())
options(stringsAsFactors = F)

wkf = "rmacro"
param = "rmacro.R"
source("/Users/Anna/Documents/workflows/rmacroRDM/R/functions.R")

## ----master-configuration, eval=T----------------------------------------
init_db(data.folder = "/Users/Anna/Google Drive/bird trait networks/",
        script.folder = "~/Documents/workflows/bird_trait_networks/", 
        spp.list_src = "D0", fileEncoding = "mac")

## ----intialise-project ----
source("~/Documents/workflows/bird_trait_networks/project_ui.R")


# SETTINGS ###############################################################



## ----setup-input.folder--------------------------------------------------
setupInputFolder(input.folder)

fcodes <- ensure_fcodes(meta.vars)

file.names <- create_file.names(c("Table.csv", "Amniote_Database_Aug_2015.csv"))

if(create.data_log){
create_data_log(file.names = c("Table.csv", "Amniote_Database_Aug_2015.csv"), 
                overwrite = T)
}

if(create.data_log){
update_data_log(overwrite = T)
}

load_sys.ref(view = F)
# Load match data.....................................................................


syn.links <- read.csv("taxo/syn.links.csv", stringsAsFactors = F)


# WORKFLOW ###############################################################


# CREATE MASTER

# Create taxo.table
taxo.dat <- unique(D0[,c("species", taxo.vars)])

# Assign spp.list from species in original dataset D0
spp.list <- createSpp.list(species = taxo.dat$species, 
                           taxo.dat = taxo.dat, 
                           taxo.vars)

# extract and order D0 into master format
D0 <- longMasterFormat(data = D0, master.vars, data.ID = "D0")
  
# create master shell
master <- list(data = newMasterData(master.vars), spp.list = spp.list, metadata = metadata)
master <- updateMaster(master, data = D0, spp.list = NULL)


# match and append processed data to master


  filename <- "D1.tidy"
  
  m <- matchObj(data.ID = "D1", spp.list = spp.list, status = "unmatched",
                data = read.csv(paste(input.folder, "csv/", filename, ".csv", sep = ""),
                                stringsAsFactors=FALSE, fileEncoding = "mac"),
                sub = "spp.list", filename = filename, 
                meta = createMeta(meta.vars)) # use addMeta function to manually add metadata.
  
  m <- processDat(m, input.folder, var.omit) %>% 
    separateDatMeta() %>% 
    compileMeta(input.folder = input.folder) %>%
    checkVarMeta(master$metadata) %>%
    dataMatchPrep()

  m <- dataSppMatch(m, syn.links = syn.links, addSpp = T)
  
  # Match data set to spp.list and process
  output <- masterDataFormat(m, meta.vars, match.vars, var.vars)
  
  write.csv(output$data, file =  "csv/D1.long.csv", row.names = F, fileEncoding = "mac")
  
  spp.list <- output$spp.list
  
  
    # dir.create(paste(output.folder, "data/", sep = ""), showWarnings = F)
    # dir.create(paste(output.folder, "data/match objects/", sep = ""), showWarnings = F)
    
  save(m, file = paste(output.folder, "data/match objects/", m$data.ID, "m.RData", sep = ""))
  
  

  
  # MERGE DATASETS
  
  master <- updateMaster(master, data = output$data, spp.list = output$spp.list)
  
  
  outliers <- read.csv("metadata/outliers.csv")
  outliers$data.ID <- "D0"
  outliers$data.ID[outliers$select == "D0"] <- "D1"

  # remove outliers
  master$data <- removeData(data = master$data, outliers = outliers)
  
  
  # remove duplicates 
  # 
  # from D0 for variable "repro.age"
  master$data <- master$data[-which(master$data$var == "repro.age" & 
                 duplicated(master$data[,c("species", "var")], fromLast = T)),]
  
  # from D1 for all other variables
  master$data <- master$data[!duplicated(master$data[,c("species", "var")]),]
  
  write.csv(master$data, file =  "analytical/master.csv", row.names = F, fileEncoding = "mac")
  
  save(master, file = paste(output.folder, "data/master.RData", sep = ""))
  
  
  # Create wide dataset ################################################################################### 
  wide <- widenMaster(vars = unique(master$data$var), species = unique(master$data$species), 
                      master = master, add.taxo = T, taxo.vars)
  
  
  write.csv(wide, file =  "analytical/master wide.csv", row.names = F, fileEncoding = "mac")
  
