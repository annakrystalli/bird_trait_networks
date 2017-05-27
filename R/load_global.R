## ---- lg-global_options ----
#' global settings
options(stringsAsFactors = F)


## ---- lg-read_pkgs ----
#' load character vector of dependency package names
source(paste(script.folder, "R/pkgs.R", sep = ""))


## ---- lg-load_pkgs ----
#' install and load dependencies through pkg "pacman"
if (!require("pacman")) install.packages("pacman")
pacman::p_unload(pacman::p_loaded(), character.only = T)
pacman::p_load(pkgs, character.only = T)


## ---- lg-load_params ----
source(paste(script.folder, "params/project_ui.R", sep = ""))
## ---- lg-load_files ----
source(paste(script.folder, "R/load_files.R", sep = ""))


## ---- lg-load_helper_functions ----
#' load helper functions
source(paste(script.folder, "R/helper_functions.R", sep = ""))


