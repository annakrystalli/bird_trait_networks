rm(list = ls())

# ---- pn-setup ----
if(exists("file_setup_path")){}else{
  file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
  source(file_setup_path)}
wkf = "phylonetworker"
param = "phylonetworker.R"
source(paste0(script.folder,"project_ui.R"))

# ---- pn-get-vg_dt ----
vg <- orderby_vg_dt(vg, mgm_types)
vg_dt <- get_vg_dt(vg, mgm_types)


# ---- pn-nn ----
# numeric - numeric edges
nn_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(x %in% c("g", "p"))}),]
nn_dat <- data[,c("species", unique(unlist(nn_vg[,1:2])))]  
nn_res <- NULL
for(i in 1:dim(nn_vg)[1]){
  print(i)
  nn_res <- rbind(nn_res, pglsPhyloCor(pair = nn_vg[i, 1:2], data = nn_dat, 
                                       tree = tree, log.vars = log.vars, 
                                       pair.type = "nn", mgm_types = mgm_types))
}
nn_res <- nn_res[order(abs(nn_res$phylocor), decreasing = T),]
if(save){
  write.csv(nn_res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                          min.n, if(log){"_log"},"_nn.csv", sep = ""),
            row.names = F)
}


# ---- pn-cc.list ----
# categorical - categorical edges (returns list)
cc_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(x == "c")}),]
cc_dat <- data[,c("species", unique(unlist(cc_vg[,1:2])))]  
cc_res.list <- vector("list", nrow(cc_vg))
for(i in 1:nrow(cc_vg)){
  print(i)
  cc_res.list[[i]] <- pglsPhyloCor(pair = cc_vg[i, 1:2], data = cc_dat, 
                                   tree = tree, log.vars = NULL, pair.type = "cc", result = "row",
                                   mgm_types = mgm_types)
}
if(save){
  save(cc_res.list, file = paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                                 min.n, if(log){"_log"},"_cc.Rdata", sep = ""))
}

# ---- pn-cc-load.list ----
load(file = paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn",
                  min.n, if(log){"_log"},"_cc.Rdata", sep = ""))

# ---- pn-cc ----
# categorical - categorical edges (list to data.frame)
cc_res <- NULL
for(i in 1:length(cc_res.list)){
  cc_res <- rbind(cc_res, get_cc.row(cc_res.list[[i]], cc_vg))
  print(i)
}
if(save){
  write.csv(cc_res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                          min.n, if(log){"_log"},"_cc.csv", sep = ""),
            row.names = F)}


# ---- pn-nc ----
# numeric - categorical edges
nc_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(c("c", "g") %in% x)}),]
nc_dat <- data[,c("species", unique(unlist(nc_vg[,1:2])))]
nc_res <- NULL
for(i in c(11:dim(nc_vg)[1])){
  nc_res <- rbind(nc_res, 
                  pglsPhyloCor(pair = nc_vg[i, 1:2], data = nc_dat, 
                               tree = tree, log.vars = log.vars, 
                               pair.type = "nc", result = "row", 
                               mgm_types = mgm_types))
  print(i)
}
if(save){write.csv(nc_res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                                 min.n, if(log){"_log"},"_nc.csv", sep = ""),
                   row.names = F)}


# ---- pn-read.type-res ----
nn_res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_nn.csv", sep = ""))
nc_res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_nc.csv", sep = ""))
cc_res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_cc.csv", sep = ""))
# ---- pn-comp.res ----
res <- rbind(nn_res, nc_res, cc_res)
if(save){
  write.csv(res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                       min.n, if(log){"_log"},"_allDT.csv", sep = ""),
            row.names = F)
}


# ---- pn-read.type-res ----
res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                      min.n, if(log){"_log"},"_allDT.csv", sep = ""))
# ---- pn-sub-res ----
res <- res[abs(res$phylocor) > cutoff & !is.na(res$phylocor),]

# ---- pn-network ----
library(rnetcarto)
net.d <- res[!is.na(res$phylocor), c("var1", "var2", "phylocor")]
net.list <- as.list(net.d)
net <- netcarto(web = net.list, seed = 1)
if(save){
  write.csv(cbind(net[[1]], modularity = net[[2]]), 
            paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                  if(log){"_log"},".csv", sep = ""),
            row.names = F)
  save(net, file = paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                         if(log){"_log"},".RData", sep = ""))
}

