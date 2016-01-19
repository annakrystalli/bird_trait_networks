rm(list=ls())

source("~/Documents/workflows/bird_trait_networks/Setup.R")


# PACKAGES & FUNCTIONS ###############################################################

source(paste(script.folder, "functions.R", sep = ""))

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

dl <- list(D1 = processDat(file = "D1.tidy.csv", label = F, taxo.dat, var.omit, input.folder,
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
                            synonyms = synonyms, taxo.table = taxo.table)
  
  write.csv(output$mdat, file =  "csv/D1.long.csv", row.names = F)
  
  master <- rbind(master, output$mdat)
  spp.list <- output$spp.list
  
  

  
  # MERGE DATASETS
  
  D1 <- output$mdat
  outliers <- read.csv("r data/outliers.csv")
  
  outRowSelect <- function(x, species, var){x[x$species == species & x$var == var, c("species", "var", "value", "ref")]} 
  outRowSelect1 <- function(x, species, var){x[x$species == species & x$var == var, c("species", "var", "value", "ref")]} 
  
out.refs <- c()
    for(i in 1:dim(outliers)[1]){
    out.refs <- rbind(cbind(outRowSelect1(D0, species = outliers[i,"species"], var = outliers[i,"var"]), dataset = "D0"),
                      cbind(outRowSelect(D1, outliers[i,1], outliers[i,2]), dataset = "D1"))
    }, 
    D1 = D1, D0 = D0)
  
  cbind(D1[D1$species == outliers$species & D1$var == outliers$var, "species", "var", "value", "ref"], D1["species", "var", "value", "ref"])
  
  duplicated(rbind(D0[,c("species", "var")], D1[,c("species", "var")]))
  
