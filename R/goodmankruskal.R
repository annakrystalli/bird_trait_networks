rm(list = ls())
wkf = "goodmankruskal"
param = "goodmankruskal.R"
file_setup_path <- "file_setup.R"
source(file_setup_path)
source(paste0(script.folder, "project_ui.R"))


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





gdk_data <- g$data[,names(mgm_types[mgm_types %in% c("c","b")])]

gtau <- GKtauDataframe(gdk_data)
diag(gtau) <- NA

m.lab <- 120
m = list(
  l = m.lab,
  r = 50,
  b = m.lab,
  t = 50,
  pad = 4
)

plot_ly(z= gtau, x = colnames(gtau), y = rownames(gtau), type = plot.type) %>%
  layout(title="GKtau among categorical variables", margins = m, 
         yaxis = list(title = "", showgrid = F), 
         xaxis = list(title = "", showgrid = F))


