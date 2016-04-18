rm(list=ls())

source("~/Documents/workflows/bird_trait_networks/setup.R")


# PACKAGES & FUNCTIONS ###############################################################

source(paste(script.folder, "functions.R", sep = ""))
source('~/Documents/workflows/Sex Roles in Birds/birds/bird app/app_output_functions.R', 
       chdir = TRUE)

require(plotly)
require(knitr)
require(RColorBrewer)


# SETTINGS ###############################################################

qcmnames = c("qc", "observer", "ref", "n", "notes")
taxo.var <- c("species", "order","family", "subspp", "parent.spp")
var.var <- c("var", "value", "data")
var.omit <- c("no_sex_maturity_d", "adult_svl_cm", "male_maturity_d")

setupInputFolder(input.folder, qcmnames)

setwd(input.folder)


# FILES ###############################################################

# source("~/Documents/workflows/bird_trait_networks/process D0.R")
# source("~/Documents/workflows/bird_trait_networks/process D1.R")

D0 <- read.csv(file = "csv/D0.csv" ,fileEncoding = "mac")                    

taxo.table <- unique(D0[, c("species", "order", "family")])

# Load match data.....................................................................
synonyms  <- read.csv("r data/synonyms.csv", stringsAsFactors=FALSE)




# WORKFLOW ###############################################################

dl <- list(D1 = processDat(file = "D1.tidy.csv", dat = NULL, label = F, taxo.dat, var.omit, input.folder,
                           observer = NULL, qc = NULL, ref = NULL, n = NULL, notes = NULL,
                           master.vname = "master.vname"))

# CREATE MASTER

# Assign spp.list from species in original dataset D3
spp.list <- data.frame(species = unique(D0$species))

# create master shell
master <- c()

# match and append processed data to master
data.ID <- "D1"

  # Create match object
  m <-  matchObj(data.ID, spp.list, data = dl[[data.ID]]$data, status = "unmatched", 
                 sub = "spp.list",
                 qcref = dl[[data.ID]]$qcref) 
  
  # Match data set to spp.list and process
  output <- matchMSToMaster(m, taxo.var = taxo.var, var.omit = var.omit, input.folder = input.folder, 
                            output.folder = output.folder, ignore.unmatched = T,
                            synonyms = synonyms, taxo.table = taxo.table,
                            trim.dat = T, retain.dup = F)
  
  write.csv(output$mdat, file =  "csv/D1.long.csv", row.names = F, fileEncoding = "mac")
  
  spp.list <- output$spp.list
  
  

  
  # MERGE DATASETS
  
  # Transform D0
  
  master <- data.frame(D0[,c("species", "order", "family")], subspp = F, parent.spp = NA,  
                D0[,c("var", "value")], data = "D0", synonyms = D0$species, 
                data.status = "original", "qc" = NA, observer = NA, ref = D0$ref, n = NA)
  
  
  D1 <- output$mdat
  outliers <- read.csv("r data/R03/outliers.csv")
  
  
  # select outliers
  
  for(i in 1:dim(outliers)[1]){
    if(outliers$select[i] == "D0"){
      D1[!(D1$species == outliers$species[i] & D1$var == outliers$var[i]),]
    }else{
      master[!(master$species == outliers$species[i] & master$var == outliers$var[i]),]
    }
  }
  
  master <- rbind(master, D1)
  
  master <- master[-which(master$var == "repro.age" & 
                 duplicated(master[,c("species", "var")], fromLast = T)),]
  
  master <- master[!duplicated(master[,c("species", "var")]),]
  
  write.csv(master, file =  "csv/master.csv", row.names = F, fileEncoding = "mac")
  
  # Create wide dataset ################################################################################### 
  wide <- widenMaster(vars = unique(master$var), species = unique(master$species), 
                      master = master, metadata = metadata)
  
  
  write.csv(wide, file =  "csv/master wide.csv", row.names = F, fileEncoding = "mac")
  
