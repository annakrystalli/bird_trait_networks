rm(list=ls())


### SETUP ##############################################################

source("~/Documents/workflows/bird_trait_networks/setup.R")
setwd(input.folder)

### LOAD RAW DATA ##############################################################

# MAY PRODUCE ERROR
# might need to check fileEncoding on your system. Try "latin1", instead of 
# "mac" to run but some characters in ref may still not display correctly 

D0 <- read.csv("csv/Table.csv", stringsAsFactors = F, fileEncoding = "mac") %>% 
  apply(2, FUN = trimws) %>% data.frame(stringsAsFactors = F)

metavar <- read.csv("csv/metavar.vnames.csv", stringsAsFactors = F, fileEncoding = "mac")

### PROCESS ##############################################################

# remove empty columns
D0 <- D0[,names(D0) %in% metavar$Table]

# rename D0 names
names(D0) <- metavar$D0[match(names(D0), metavar$Table)]

# correct error in var entry for Branta sandvicensis
D0[D0$species == "Branta sandvicensis" & D0$var == "sandvicensis", "var"] <- "Age at first reproduction"

### assign code variable name to D0$var ###
D0$var <- vnames$code[match(D0$var, vnames$D0)]

# hyphenate species names
D0$species <- gsub(" ", "_", D0$species)
D0$value[D0$var == "repro.age"] <- as.numeric(D0$value[D0$var == "repro.age"]) * 365

# write processed data
write.csv(D0, file = "csv/D0.csv", row.names = F, fileEncoding = "mac")                    
write.csv(sort(unique(D0$var)), file = "r data/var.vnames.csv", row.names = F, fileEncoding = "mac")                    
                 