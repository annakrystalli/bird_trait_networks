rm(list=ls())


### SETUP ##############################################################

source("~/Documents/workflows/bird_trait_networks/setup.R")

qcmnames = c("qc", "observer", "ref", "n", "notes")
taxo.var <- c("species", "order","family", "subspp", "parent.spp")
var.var <- c("var", "value", "data")
var.omit <- c("no_sex_maturity_d", "adult_svl_cm", "male_maturity_d")

# prepare match folders
setupInputFolder(input.folder, qcmnames)
setwd(input.folder)

# FILES ###############################################################


D1 <- read.csv(paste("csv/", "Amniote_Database_Aug_2015.csv", sep = ""), stringsAsFactors = F, 
               fileEncoding = "mac") %>% 
  apply(2, FUN = trimws) %>% data.frame(stringsAsFactors = F)
                
 
# WORKFLOW ###############################################################



##### Initial D1 processing ###########


### subset D1 by class to aves and process
D1 <- D1[D1$class == "Aves",]
D1 <- D1[,!apply(D1, 2, FUN = function(x){all(is.na(x))})]
D1$species <- paste(D1$genus, D1$species, sep = "_")

# convert adult_svl_cm from cm to mm
D1$adult_svl_cm <- as.numeric(D1$adult_svl_cm) * 10

# separate and save taxonomic data
taxo <- D1[, names(D1) %in% c("class", "order", "family", "genus","species" ,"common_name")]
D1 <- D1[, !names(D1) %in% c("class", "order", "family", "genus","common_name")]
write.csv(taxo, file = "r data/taxo.csv", row.names = F)

# write proceesed D1 file, explored in R01
write.csv(D1, file = "csv/D1.csv", row.names = F, fileEncoding = "mac")

D1.names <- names(D1)

# Post R01 D1 data ###########################################################################

D1 <- read.csv(file = paste("csv/","D1.csv", sep = ""),
               fileEncoding = "mac")

# correct typo in "male_maturity_d"
D1[D1$species == "Rhea_americana","male_maturity_d"] <- D1[D1$species == "Rhea_americana","male_maturity_d"]/10
D1$repro.age.diff <- D1$female_maturity_d - D1$male_maturity_d
D1 <- D1[, !names(D1) %in% var.omit]

D1 <- codeVars(dat = D1, data.ID = "D1", metadata = metadata, vnames = vnames)

### SAVE ###
write.csv(D1, paste(getwd(),"/csv/D1.tidy.csv", sep = ""), row.names = F)
#_______________________________________________________________


##### Process D1 refs #####
D1.ref <- read.csv(file = paste("ref/","D1.csv", sep = ""),
                   fileEncoding = "mac")
D1.ref <- D1.ref[D1.ref$species %in% D1$species,]

# combine references for derived repro.age.diff
D1.ref$repro.age.diff <- paste(D1.ref$female_maturity_d, D1.ref$male_maturity_d, sep = "; ")
D1.ref$repro.age.diff[is.na(D1.ref$female_maturity_d) | is.na(D1.ref$male_maturity_d)] <- NA
D1.ref <- D1.ref[, !names(D1.ref) %in% var.omit]

# code variables
D1.ref <- codeVars(dat = D1.ref, data.ID = "D1", metadata = metadata, vnames = vnames)


### SAVE ###
write.csv(D1.ref, paste(getwd(), "/ref/D1.tidy.csv", sep = ""), row.names = F)
#_________________________________________________________________



##### Process D1 n #####
D1.n <- read.csv(file = paste(getwd(),"/n/","D1.csv", sep = ""),
                 fileEncoding = "mac")
D1.n <- D1.n[D1.n$species %in% D1$species,]

# remove range variables
D1.n <- D1.n[,-c(grep("min", names(D1.n)), grep("max", names(D1.n)))]
# correct var names
names(D1.n) <- gsub("count_", "", names(D1.n))
names(D1.n) <- unlist(lapply(names(D1.n), FUN = function(x, pattern){
  m <- grep(pattern = x, x = pattern, value = T)
  if(length(m) == 0){m <- x}
  if(length(m) == 2){m <- m[2]}
  m}, 
pattern = D1.names))

# combine references for derived repro.age.diff
D1.n$repro.age.diff <- D1.n$female_maturity + D1.n$male_maturity
D1.n[D1.n == 0] <- NA
D1.n <- D1.n[, !names(D1.n) %in% var.omit]

# code variables
D1.n <- codeVars(dat = D1.n, data.ID = "D1", metadata = metadata, vnames = vnames)

### SAVE ###
write.csv(D1.n, paste(getwd(),"/n/D1.tidy.csv", sep = ""), row.names = F)
#_________________________________________________________________

