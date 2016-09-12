
# ---- full.vg ----
g <- m1DataPrep(g, phylo.match)
g$data <- g$data[, intersect(names(g$data),
                             names(g$mgm_type)[g$mgm_type != "p"])]
g <-  aggDat(g, min.cat)
g <- rm_sing(g)

#' vector of no. of levels for categorical/binomial variables. 1 for all 
#' other data types
#+ lev.mod
g$lev <- apply(g$data, 2 , FUN= function(x){length(unique(na.omit(x)))})
g$lev[g$mgm_types[names(g$lev)] %in% c("g", "p")]  <- 1 

if(log){
g <- logG(g)}
