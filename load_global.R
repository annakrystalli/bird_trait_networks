#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---

#' global settings

# ---- load_global ----
options(stringsAsFactors = F)

# dependencies
source(paste(script.folder, "pkgs.R", sep = ""))

# install dependencies
if(install.pkgs){
  source(paste(script.folder, "install_dependencies.R", sep = ""))}

# load dependencies
source(paste(script.folder, "load_dependencies.R", sep = ""))

# load files
source(paste(script.folder, "load_files.R", sep = ""))


