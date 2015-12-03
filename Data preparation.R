rm(list=ls())

options(stringsAsFactors = F)

### SETTINGS ##############################################################

codewd <- "~/Documents/worflows/bird_trait_networks/"
datawd <- "~/"

output.folder <- "/Users/Anna/Google Drive/bird trait networks/outputs/"
input.folder <- "/Users/Anna/Google Drive/bird trait networks/inputs/data/"
setwd(paste(input.folder, "csv", sep = ""))





### LOAD DATA ##############################################################

D0 <- read.csv("Tabla.csv", stringsAsFactors = F, fileEncoding = "mac") 
# might need to check fileEncoding on your system. Try "latin1", instead of 
# "mac" to run but some characters in ref might not display correctly 

