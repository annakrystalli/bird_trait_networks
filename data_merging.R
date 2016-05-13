# SETUP ###############################################################


rm(list=ls())

source("~/Documents/workflows/bird_trait_networks/Setup.R")
setwd(input.folder)

# PACKAGES & FUNCTIONS ###############################################################

source(paste(script.folder, "functions.R", sep = ""))


# install.packages("prodlim")
require(prodlim)

# SETTINGS ###############################################################





# Determine outliers

D0 <- read.csv(file = paste(input.folder, "csv/","D0.csv", sep = "")
               ,fileEncoding = "mac")                    
D1 <- read.csv(file = paste(input.folder, "csv/","D1.long.csv", sep = ""),
               fileEncoding = "mac") 

### Remove duplicates from D0. 
D0 <- D0[!duplicated(D0[,c("species", "var")]),]

## Check 

D <- rbind(cbind(D0[,c("species", "var", "value")], dataset = "D0"),
           cbind(D1[,c("species", "var", "value")], dataset = "D1"))



# Deal with outliers
outliers.done <- read.csv("~/Google Drive/bird trait networks/inputs/data/r data/outliers done.csv", stringsAsFactors=FALSE)






dups <- D[duplicated(D[,c("species", "var")]),c("species", "var")]

dat <- D[row.match(dups, D[,c("species", "var")]),]
  
dat <- cbind(dups, D0.value = D0[row.match(dups, D0[,c("species", "var")]), "value"],
      D1.value = D1[row.match(dups, D1[,c("species", "var")]), "value"])

D0[D0$species == dups$species & D0$var == dups$var, ]


if(D0var != "bill.len"){
  df <- cbind(dat$D0.value, dat$D1.value)
  outliers <- apply(df, 1, FUN = 
                      function(x){any(abs(diff(x))/x[which(x == min(x))] > 2)})
  outlie.df <- dat[!outliers,]
  
  outlie.df <- rbind(outlie.df, cbind(species = spp, var = D0var, df)[outliers,])
} 
