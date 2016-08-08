#' ---
#' output: 
#' html_document:
#' theme: paper
#' ---

#' # Project Interface
#' ***
#' <br>

#+ md-setup, echo=F, message=F, warning=F
rm(list = ls())
require(knitr)
opts_chunk$set(echo=TRUE, message=F, warning=F, cache=FALSE)
#'
#+ file-setup, echo=F
file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
#'
#+ read_chunks, echo = F
read_chunk(paste(script.folder, "file_setup.R", sep = ""))
read_chunk(paste(script.folder, "load_dependencies.R", sep = ""))
read_chunk(paste(script.folder, "load_files.R", sep = ""))
read_chunk(paste(script.folder, "load_environment.R", sep = ""))
#'

#' ## **get files**
#' Before starting download the data and github files from osf:
#' [`bird traits`](https://osf.io/h7kgy/). You'll also need to download the [`rmacroRDM`](https://github.com/annakrystalli/rmacroRDM) folder.
#' #' You can also access authorised data directly through r by [installing googledrive on your computer](https://tools.google.com/dlpage/drive). You will be then be able to specify direct paths to the data in the google drive folder.

#'
#'***
#'
#' <br>
#' 
#' ## **setup `file system`**
#' 
#' ### - edit the information in `file_setup.R`
#' Link the workflow to your local copies of the folders comprising the **project file system**:
#' 
#' #### current file setup
#' ##### **data input folders:** - link to **google drive folder**
#+ data-setup
#' ##### **script folders:** - link to cloned **`github repo folders`**
#+ script-setup

#' 
#' <br>
#' 
#' ### - source `file_setup.R`.
#' full path to the file is required
#+ file-setup, eval=F

#' **check objects in the working environment:** 
#' *`r ls()`*
#' 

#' 
#' ***
#' 
#' <br>
#' 
#' ## **initialise project environment**
#' #### - require package installation?
#+ install-pgks?
install.pkgs <- F
#'
#'
#' ### - source **`load_global.R`**
#+ init_global, echo = T
source(paste(script.folder, "load_global.R", sep = ""))


#'<br>
#'
#' #### > source contents
#' 
#' ##### ***load dependencies `load_dependencies.R`***
#+ load_dependence, results = 'asis', eval= F
#'
#' **packages loaded:** *`r pkgs`*
#' 
#'
#' ##### ***load global files `load_files.R`***
#' load global files required across project
#+ require-dplyr, results = 'asis', eval= F
#'

#+ load_global_files, results = 'asis', eval= F
#'
#' **check global setup:** 
#' *`r ls()`*
#' 
#' ***
#' <br>
#'
#' ## **initialise analysis environment**

#' ### - set analyses parameters
#' Set analyses parameters by specifying the workflow (`wkf`) and `param` file name in `script_folder, "/params"`
#+ analyses_params
wkf = "mgm"
param = "phylonetworker.R"
#'  
#'<br>
#'
#' ### - source analyses environment
#+ source-analyses
source(paste(script.folder, "load_environment.R", sep =""))
read_chunk(paste(script.folder,"params/", param, sep = ""))
#'

#'<br>
#'
#' ####     > source contents 
#'
#' ##### ***load analysis parameters: specified `params/` `.R` file***
#' full file name (includind `.R`) required.
#+ source-params, results = 'asis'
#'

#' ***current `/params/... .R` file: ***  `r param` 
#+ params, results = 'asis'
#'

#'<br>
#'
#' ##### ***load analysis environment***
#+ load_phylonetworker, results = 'asis'
#'

