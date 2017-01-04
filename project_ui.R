
# ----- pui-init_global ----
source(paste(script.folder, "load_global.R", sep = ""))

# ----- pui-source-environment ----
source(paste(script.folder, "load_environment.R", sep =""))
read_chunk(paste(script.folder,"params/", param, sep = ""))
