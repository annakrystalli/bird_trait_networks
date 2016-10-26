#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---

#' global settings

## ---- load_global ----
options(stringsAsFactors = F)

# dependencies
source(paste(script.folder, "pkgs.R", sep = ""))

# install and load dependencies through pkg "pacman"
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pkgs, character.only = T)

# load files
source(paste(script.folder, "load_files.R", sep = ""))

# load helper functions
source(paste(script.folder, "helper_functions.R", sep = ""))


