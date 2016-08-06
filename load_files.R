#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---
#' 
#' ### load metadata
# ---- require-dplyr ----
require(dplyr)

# ---- load_global_files
metadata <- read.csv(paste(input.folder, "metadata/","metadata.csv", sep = ""), 
                     stringsAsFactors = F, fileEncoding = "mac") %>% 
  apply(2, FUN = trimws) %>% data.frame(stringsAsFactors = F)

vnames <- read.csv(paste(input.folder, "metadata/","vnames.csv", sep = ""), 
                   stringsAsFactors = F, fileEncoding = "mac") %>% 
  apply(2, FUN = trimws) %>% data.frame(stringsAsFactors = F)
vnames[vnames == ""] <- NA


