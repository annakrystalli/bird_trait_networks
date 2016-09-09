# ----- req-knit ----
require(knitr)

# ----- read_chunks -----
read_chunk(paste(script.folder, "file_setup.R", sep = ""))
read_chunk(paste(script.folder, "pkgs.R", sep = ""))
read_chunk(paste(script.folder, "load_dependencies.R", sep = ""))
read_chunk(paste(script.folder, "load_global.R", sep = ""))
read_chunk(paste(script.folder, "helper_functions.R", sep = ""))
read_chunk(paste(script.folder, "load_files.R", sep = ""))
read_chunk(paste(script.folder, "load_environment.R", sep = ""))

# ----- install-pgks? ----
install.pkgs <- F

# ----- init_global ----
source(paste(script.folder, "load_global.R", sep = ""))

# ----- source-analyses
source(paste(script.folder, "load_environment.R", sep =""))
read_chunk(paste(script.folder,"params/", param, sep = ""))
