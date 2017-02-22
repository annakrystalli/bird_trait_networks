## ----global-setup, echo = F----------------------------------------------
rm(list=ls())
options(stringsAsFactors = F)

wkf = "rmacro"
param = "rmacro.R"

if(exists("file_setup_path")){}else{
  file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
  source(file_setup_path)}
source(paste0(script.folder,"project_ui.R"))

launch_sr_configurator(sr_configurator, file_setup_path)

## ----master-configuration, eval=T----------------------------------------
init_db(data.folder = data.folder,
        script.folder = getwd(), 
        spp.list_src = "D0", fileEncoding = "UTF-8")

# SETTINGS ###############################################################
# ---- load-sy.ref ----
load_sys.ref(view = F)

# WORKFLOW ###############################################################

process_file.system(save.taxo = T)

# Load match data.....................................................................
syn_links <- read.csv(paste0(input.folder, "taxo/syn_links.csv"), stringsAsFactors = F)

# Assign spp.list from species in original dataset D0
spp.list <- createSpp.list(syn_links = syn_links)
  
# create master shell
master <- create_master(file.name = get_file.names()["D0"], 
                        spp.list = spp.list)
## ----create-m------------------------------------------------------------

filename <- get_file.names()["D1"]

m <- matchObj(file.name = filename,
              spp.list = master$spp.list,
              sub = "spp.list") # use addMeta function to manually add metadata.


## ----process-m-----------------------------------------------------------
m <- m %>% 
  separateDatMeta() %>% 
  compileMeta() %>%
  checkVarMeta() %>%
  dataMatchPrep()

## ----data-spp-match------------------------------------------------------
m <- dataSppMatch(m, addSpp = T)
save(m, file = paste(output.folder, "data/match objects/", m$dcode, "m.rda", sep = ""))

## ----output--------------------------------------------------------------
output <- masterDataFormat(m)

## ----merge-to-master-----------------------------------------------------
master <- updateMaster(master, output = output)

  
  spp.list <- output$spp.list
  
  # ---- post-process ----
  # remove outliers
  outliers <- read.csv(paste0(ds$input.folder, "metadata/outliers.csv"))
  outliers$data.ID <- "D0"
  outliers$data.ID[outliers$select == "D0"] <- "D1"
  master$data <- removeData(data = master$data, outliers = outliers)
  
  # remove duplicates 
  # from D0 for variable "repro.age"
  master$data <- master$data[-which(master$data$var == "repro.age" & 
                 duplicated(master$data[,c("species", "var")], fromLast = T)),]
  
  # from D1 for all other variables
  master$data <- master$data[!duplicated(master$data[,c("species", "var")]),]
  
  # ---- save-outputs ----
  write.csv(master$data, file =  paste0(input.folder, "analytical/master.csv"), 
            row.names = F, fileEncoding = ds$fileEncoding)
  save(master, file = paste(output.folder, "data/master.rda", sep = ""))
  
  
  # ---- create wide dataset ----
  wide <- widenMaster(master, vars = unique(master$data$var), species = unique(master$data$species), 
                      add.taxo = T)
  
  
  write.csv(wide, file =   paste0(input.folder, "analytical/master wide.csv"), row.names = F, 
            fileEncoding = ds$fileEncoding)
  
