# ---- init-an ----
source(paste(script.folder, "project_ui.R", sep = ""))


num <- unique(unlist(num.vg[,1:2]))
bin <- unique(unlist(bin.vg[,1:2]))
int <- unique(unlist(int.vg[,1:2]))
cats <- unique(unlist(cat.vg[,1:2]))


num.id <- 1:15
int.id <- 1:5
bin.id <- 1:10
cat.id <- 6:10

vars <- c(num[num.id], int[int.id], bin[bin.id], cats[cat.id])

type <- c(rep("g", length(num.id)),
          rep("p", length(int.id)),
          rep("c", length(c(bin.id, cat.id))))

data <- m1DataPrep(data = data, datType = unique(metadata$type), phylo.match, log, log.vars)

# aggregate or remove categories with < min.cat frequency
#________________________________________________________________________

# tabulate each categorical variable
tabs <- apply(cat.dat[,!names(cat.dat) %in% c("species", "synonyms")],
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

g <- list(data = data, meta = metadata)

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

for(var in names(other)){
  
  g <- aggregateCats(g, var, agg.cats = other[[var]])
}

for(var in names(nas)){
  
  g <- naCats(g, var, agg.cats = nas[[var]])
}


# Select data and fit model ##########################################


dat <- g$data[,c(vars)]


lev.mod <- apply(dat, 2 , FUN= function(x){length(unique(na.omit(x)))})
lev.mod[type %in% c("g", "p")] <- 1 




test <- mgmfit(dat, type, lev.mod, lambda.sel = "EBIC", rule.reg = "OR", method = "glm",
               missings = "casewise.zw")

# define variable types
groups_typeV <- list("Gaussian"=which(type=='g'), 
                     "Categorical"=which(type=='c'))

# pick some nice colors
group_col <- c("#72CF53", "#ED3939")


qgraph(test$wadj, 
       vsize=3.5, 
       esize=5, 
       edge.color = test$edgecolor, 
       color=group_col,
       border.width=1.5,
       border.color="black",
       groups=groups_typeV,
       nodeNames=names(dat),
       legend=TRUE, 
       legend.mode="style2",
       legend.cex=0.3)


