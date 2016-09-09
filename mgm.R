#' ---
#' output: 
#' html_document:
#' theme: paper
#' self_contained: false
#' ---
#' 
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

vars <- ms_vars[1:(as.integer(length(ms_vars)*v.p))]
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

source(paste(script.folder, "mgm_dataprep.R", sep = ""))

#' impute-mgm
source(paste(script.folder, "impute_mgm.R", sep = ""))


#' proportion of data imputed: `r `
sum(is.na(g$data))/prod(dim(g$data))

#' ### fit mgm
#+ fit-mgm, cache = T
mgm_mod <- mgmfit(g$imp_data, type = g$mgm_types[names(g$data)], 
               lev = g$lev[names(g$data)], lambda.sel = "EBIC", 
               rule.reg = "OR", method = "glm", missings = "casewise.zw")

#' ### fit netcarto
#+ fit-netcarto, cache = T
library(rnetcarto)
net <- rnetcarto::netcarto(mgm_mod$wadj)


save(g, mgm_mod, net, file = paste(output.folder, "data/mgm/", an.ID, 
                                   "-", v.p, "-", s.p, "-", min.n, 
                                   "-", min.cat, ".Rdata", sep = ""))

#' ## plot
#' define variable types
#+ group-type
groups_typeV <- split(as.numeric(net[[1]]$name), f = factor(net[[1]]$module))

#' pick some nice colors
#+ group-col
group_col <- rainbow(length(groups_typeV))

#' ### plot
#+ plot, fig.width = 12, fig.height = 12
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
       legend.cex=0.1,
       details = T)


qgraph(test$wadj, 
       vsize=2,
       layout="spring",
       esize=5,
       edge.color = test$edgecolor, 
       color=group_col,
       border.width=1.5,
       border.color="black",
       minimum = "sig",
       bonf = T,
       sampleSize = prod(dim(imp.dat)),
       groups=groups_typeV,
       nodeNames=names(dat),
       overlay = T,
       legend=TRUE, 
       legend.mode="style2",
       legend.cex=0.112,
       details = T)





