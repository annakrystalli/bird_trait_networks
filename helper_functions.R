#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---

## ---- helper_functions ----
numerise <- function(x){if(all(grepl('^[0-9.]+$', x))) as.numeric(x) else x}

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# matches variable names to codes
codeVars <- function(dat, data.ID, metadata = metadata, vnames = vnames){
  # code new variables added to dataset using metadata table
  names(dat)[match(metadata$orig.vname[which(metadata$orig.vname %in% names(dat))], names(dat))] <- 
    metadata$code[which(metadata$orig.vname %in% names(dat))]
  
  # code variable names which match directly to variables in master using vnames tables
  names(dat)[match(vnames[,data.ID][which(vnames[,data.ID] %in% names(dat))], names(dat))] <- 
    vnames$code[which(vnames[,data.ID] %in% names(dat))]  
  
  if(any(is.na(names(dat)))){stop("error in coding variable. No matching code found")}
  
  return(dat)}


# function to aggregate levels below min.n into others
aggregateCats <- function(g, var, agg.cats){
  
  scores  <- g$meta[g$meta$code == var, "scores"]
  levels  <- g$meta[g$meta$code == var, "levels"]
  
  add <- max(as.numeric(strsplit(g$meta[g$meta$code == var, "scores"],
                                 ";")[[1]])) + 1
  
  g$meta[g$meta$code == var, "scores"] <- paste(scores, ";", add, sep = "")
  g$meta[g$meta$code == var, "levels"] <- paste(levels, ";other", sep = "")
  
  g$data[,var][g$data[,var] %in% agg.cats] <- add
  
  return(g)
  
}

# function to convert levels to NA
naCats <- function(g, var, agg.cats){
  
  g$data[,var][g$data[,var] %in% agg.cats] <- NA
  
  return(g)
  
}
