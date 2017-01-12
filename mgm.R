
#' # mixed graphical model analysis of bird trait data.
#' ***
#' <br>
#'
#' ## initialise project
# ---- init-an ----
rm(list = ls())
wkf = "mgm"
param = "mgm.R"
file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
source("~/Documents/workflows/bird_trait_networks/project_ui.R")

#+ define-dirs, include = F
knitr::opts_chunk$set(echo = T, warning = F, message = F)

#' ## process data to `g` object
if(is.null(dist.types)){
  vars <- ms_vars
  }else{
  vars <- ms_vars[mgm_types %in% dist.types]
  }
vars <- vars[1:(as.integer(length(vars)*v.p))]
spps <- spp.ranks$spp[1:(as.integer(nrow(spp.ranks)*s.p))]
spp <- data$species %in% spps
rownames(data) <- data$species
g <- list(data = data[spp, vars], 
          meta = metadata, 
          mgm_types = mgm_types, meta_types = meta_types, 
          log.vars = log.vars, 
          spp.list = data.frame(species = data[spp, c("species")]),
          v.p = v.p, s.p = s.p, phylo.match = phylo.match[spp,],
          tree = drop.tip(tree, setdiff(tree$tip.label, data$species[spp])))

#' ## prepare data for imputation
source(paste(script.folder, "mgm_dataprep.R", sep = ""))

#' ## impute data
source(paste(script.folder, "impute_mgm.R", sep = ""))


#' proportion of data imputed: `r sum(is.na(g$data))/prod(dim(g$data))`


#' ## fit mgm
#' fit mgm on imputed data
#+ fit-mgm, cache = T, results = "hide"
mgm_mod <- mgmfit(g$imp_data, type = g$mgm_types[names(g$data)], 
               lev = g$lev[names(g$data)], lambda.sel = "EBIC", 
               rule.reg = "OR", method = "glm", missings = "casewise.zw")

#' ## fit netcarto
#' fit network on mgm adjacency matrix
#+ fit-netcarto, cache = T
library(rnetcarto)
net <- rnetcarto::netcarto(mgm_mod$wadj)

#' #### save analysis objects for supplying to **mgm_viz.Rmd**
save(g, mgm_mod, net, file = paste(output.folder, "data/mgm/", an.ID, 
                                   "-", v.p, "-", s.p, "-", min.n, 
                                   "-", min.cat, ".Rdata", sep = ""))

#' ## plot network
#' define modules
#+ group-type
groups_typeV <- split(as.numeric(net[[1]]$name), f = factor(net[[1]]$module))

#' pick some nice colors
#+ group-col
group_col <- rainbow(length(groups_typeV))

#' ### plot
#' example of network plot using `qgraph`
#+ plot, fig.width = 12, fig.height = 12
require(qgraph)
qgraph(mgm_mod$wadj, 
       vsize=2,
       layout="spring",
       esize=5,
       edge.color = mgm_mod$edgecolor, 
       color=group_col,
       border.width=1.5,
       border.color="black",
       groups=groups_typeV,
       nodeNames=names(g$imp_data),
       overlay = T,
       legend=TRUE, 
       legend.mode="style2",
       legend.cex=0.4,
       details = T)

#' or use [**tools/mgm_viz.Rmd**](https://github.com/annakrystalli/bird_trait_networks/blob/master/reports/results/mgm_viz.Rmd). 
#' See [**example**](http://annakrystalli.github.io/bird_trait_networks/mgm_viz.nb.html)

