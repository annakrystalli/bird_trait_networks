rm(list=ls())

source("~/Documents/workflows/bird_trait_networks/setup.R")


# PACKAGES & FUNCTIONS ###############################################################

# source rmacroRDM functions
source(paste(script.folder, "functions.R", sep = ""))
source(paste(script.folder, "wideData_function.R", sep = ""))

require(plotly)
require(knitr)
require(RColorBrewer)


# SETTINGS ###############################################################
# master settings
var.vars <- c("var", "value", "data.ID")
match.vars <- c("synonyms", "data.status")
meta.vars = c("qc", "observer", "ref", "n", "notes")
master.vars <- c("species", match.vars, var.vars, meta.vars)

# spp.list settings
taxo.vars <- c("genus", "family", "order")

# var.omit settings
var.omit <- c("no_sex_maturity_d", "adult_svl_cm", "male_maturity_d")

setupInputFolder(input.folder, meta.vars)

setwd(input.folder)


# FILES ###############################################################

# source("~/Documents/workflows/bird_trait_networks/process D0.R")
# source("~/Documents/workflows/bird_trait_networks/process D1.R")

D0 <- read.csv(file = "csv/D0.csv" ,fileEncoding = "mac")                    


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
  
  
  outliers <- read.csv("r data/R03/outliers.csv")
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
  
  write.csv(master$data, file =  "csv/master.csv", row.names = F, fileEncoding = "mac")
  
  save(master, file = paste(output.folder, "data/master.RData", sep = ""))
  
  
  # Create wide dataset ################################################################################### 
  wide <- widenMaster(vars = unique(master$data$var), species = unique(master$data$species), 
                      master = master, add.taxo = T, taxo.vars)
  
  
  write.csv(wide, file =  "csv/master wide.csv", row.names = F, fileEncoding = "mac")
  
