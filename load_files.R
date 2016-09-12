#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---
#' 
#' ### load metadata
# ---- require-dplyr ----
metadata <- read.csv(paste(input.folder, "metadata/","metadata.csv", sep = ""), 
                     stringsAsFactors = F, fileEncoding = "mac")
metadata <- data.frame(apply(metadata, MARGIN = 2, FUN = trimws), stringsAsFactors = F)

vnames <- read.csv(paste(input.folder, "metadata/","vnames.csv", sep = ""), 
                   stringsAsFactors = F, fileEncoding = "mac") 
vnames <- data.frame(apply(vnames, 2, FUN = trimws), stringsAsFactors = F)  
vnames[vnames == ""] <- NA

