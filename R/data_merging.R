
## ----global-setup, echo = F----------------------------------------------
rm(list=ls())
options(stringsAsFactors = F)

wkf = "rmacro"
param = "rmacro.R"
file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
source("~/Documents/workflows/bird_trait_networks/project_ui.R")

setwd(input.folder)



# SETTINGS ###############################################################

## ----master-configuration, eval=T----------------------------------------
init_db(spp.list_src = "D0")

## ----setup-input.folder--------------------------------------------------
setupInputFolder(input.folder)

fcodes <- ensure_fcodes(meta.vars)

file.names <- create_file.names(c("Table.csv", "Amniote_Database_Aug_2015.csv"))


load_sys.ref(view = T)


## ----load-syn.links------------------------------------------------------
syn.links <- read.csv(text=getURL("https://raw.githubusercontent.com/annakrystalli/rmacroRDM/master/data/input/taxo/syn.links.csv", 
                                  ssl.verifypeer = FALSE), header=T)

## ----process-csvs, warning=FALSE-----------------------------------------

process_file.system(file.names, fcodes)


## ----create-spp.list-----------------------------------------------------
spp.list <- createSpp.list(species = NULL, 
                           taxo.dat = NULL, 
                           spp.list_src = spp.list_src)


## ----create-master-------------------------------------------------------
master <- create_master(spp.list)

## ----create-m------------------------------------------------------------

filename <- file.names[file.names == "Amniote_Database_Aug_2015.csv"]

m <- matchObj(file.name = filename,
              spp.list = master$spp.list,
              sub = "spp.list") # use addMeta function to manually add metadata.


## ----process-m-----------------------------------------------------------
m <- m %>% 
  separateDatMeta() %>% 
  compileMeta(input.folder = input.folder) %>%
  checkVarMeta(master$metadata) %>%
  dataMatchPrep()

## ----data-spp-match------------------------------------------------------
m <- dataSppMatch(m, syn.links = syn.links, addSpp = T)

## ----output--------------------------------------------------------------
output <- masterDataFormat(m, meta.vars, match.vars, var.vars)

## ----merge-to-master-----------------------------------------------------
master <- updateMaster(master, output = output)








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
