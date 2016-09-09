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
  if(length(unique(g$data[,var])) == 1){
    print(paste("var", var, "collapsed to single category - removed"))
    g$data <- g$data[, names(g$data) != var]
  }
  return(g)
  
}

# function to convert levels to NA
naCats <- function(g, var, na.cats){
  
  g$data[,var][g$data[,var] %in% na.cats] <- NA

  if(all(is.na(g$data[,var]))){
    print(paste("vars", var, "all NA - removed"))
  g$data <- g$data[, names(g$data) != var]
  g$all.na <- c(g$all.na, var)}
  
  return(g)
  
}


# aggregate or remove categories with < min.cat frequency
# tabulate each categorical variable
aggDat <- function(g, min.cat = 8) {
  tabs <- apply(g$data[,g$mgm_types[names(g$data)] %in% "c"],
                2, FUN = table) 
  # identify catecorical variables that need aggregating
  agg.var <- sapply(tabs, FUN = function(x){any(x < min.cat)})
  # identify levels below minimum n
  levs <- tabs %>% lapply(FUN = function(x){names(which(x < min.cat))})
  # identify variables in which aggregating small n categories would yield enough
  #  data for 'other' category above minimum n
  agg <- mapply(FUN = function(tabs, levs){sum(tabs[levs])} > min.cat, tabs, levs)
  
  # create indicator vector to apply transformations
  other <- levs[agg & agg.var]
  nas <- levs[(!agg) & agg.var]
 
  for(var in names(other)){
    g <- aggregateCats(g, var, agg.cats = other[[var]])
  }
  for(var in names(nas)){
    g <- naCats(g, var, na.cats = nas[[var]])
  }
  return(g)
}


rm_sing <- function(g) {
  cats <- intersect(names(g$mgm_types)[g$mgm_types == "c"], 
                    names(g$data))
  
  rm <- cats[lapply(cats, function(x, data){table(data[,x])}, 
                       data = g$data) %>%
                  sapply(function(x)length(x) == 1)]
  
  if(length(rm) == 0){return(g)}else{
    g$data <- g$data[, !names(g$data) %in% rm]
    return(g)
  }
}