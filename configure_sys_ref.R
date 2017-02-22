
setwd(script.folder)

create_data_log(fcodes, file.names = c("Table.csv", "Amniote_Database_Aug_2015.csv"), 
                overwrite = F) 
data_log <- update_data_log(overwrite = F)

vnames <- create_vnames(overwrite = F)
vnames <- update_vnames(save =T)



source('shiny/data_log.R')
