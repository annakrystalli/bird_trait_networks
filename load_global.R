#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---

#' global settings

## ---- lg-global_options ----
options(stringsAsFactors = F)

#' dependencies
# ## ---- lg-read_pkgs ----
source(paste(script.folder, "pkgs.R", sep = ""))

#' install and load dependencies through pkg "pacman"
## ---- lg-load_pkgs ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(pkgs, character.only = T)
eval(parse(text = getURL("https://raw.githubusercontent.com/annakrystalli/rmacroRDM/master/R/functions.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://raw.githubusercontent.com/annakrystalli/rmacroRDM/master/R/wideData_function.R", ssl.verifypeer = FALSE)))

#' load files
## ---- lg-load_files ----
source(paste(script.folder, "load_files.R", sep = ""))

#' load helper functions
## ---- lg-load_helper_functions ----
source(paste(script.folder, "helper_functions.R", sep = ""))


