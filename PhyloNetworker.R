# ---- pn-clear-wkf ----
.rs.restartR()
rm(list = ls(all.names = T))

# ---- pn-setup ----
if(exists("file_setup_path")){}else{
  file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
  source(file_setup_path)}
wkf = "phylonetworker"
param = "phylonetworker.R"
source(paste0(script.folder,"project_ui.R"))

# ---- pn-get-TD ----
if(!file.exists(paste(input.folder, "taxo/TD_", 
                      an.ID, 
                      if(sex_dm){"_sdm"}else{},
                      if(var.delete){"_vd"}else{},
                      ".Rdata", sep = ""))){
  # get TD and extract into df
  td <- getTD(data, tree, all.vars = F)
  TDdf <- extractTD(td, vars = ms_vars)
  # get phylosig
  phylosig <- mapply(FUN = function(x, species,tree){phylosig(tree, setNames(x, species),
                                                  method = "lambda")$lambda}, 
         x = select(data, -species),
         MoreArgs = list(species = data$species, tree = tree))
  TDdf$phylosig <- NA
  TDdf$phylosig[match(TDdf$var, names(phylosig))] <- phylosig
  save(td, TDdf, file = paste(input.folder, "taxo/TD_", 
                              an.ID, 
                              if(sex_dm){"_sdm"}else{},
                              if(var.delete){"_vd"}else{},
                              ".Rdata", sep = ""))}else{
                                load(file = paste(input.folder, "taxo/TD_", 
                                                  an.ID, 
                                                  if(sex_dm){"_sdm"}else{},
                                                  if(var.delete){"_vd"}else{},
                                                  ".Rdata", sep = ""))}

# ---- pn-get-vg_dt ----
vg <- orderby_vg_dt(vg, mgm_types)
vg_dt <- get_vg_dt(vg, mgm_types)

# ---- pn-alco-correct ----
if(file.exists(paste(output.folder, "data/phylocors/", "alco_list", 
                     min.n, if(log){"_log"},"_.Rdata", sep = ""))){
  load(paste(output.folder, "data/phylocors/", "alco_list", 
             min.n, if(log){"_log"},"_.Rdata", sep = ""))
}else{
  alco_list <- alco_correct_data(data, alco.vars, tree, mgm_types, pair.type = "nn")
  save(alco_list, file = paste(output.folder, "data/phylocors/", "alco_list", 
                               min.n, if(log){"_log"},"_.Rdata", sep = ""))
}
alco.vars <- names(alco_list$alco.complete)[alco_list$alco.complete]

# ---- pn-set-parallel ----
#cl <- makeCluster(3) 
#registerDoParallel(cl)
#clusterExport(cl, varlist = "pkgs")

# ---- pn-nn -------------------------------------------------------------------
# numeric - numeric edges
nn_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(x %in% c("g", "p"))}),]
nn_dat <- data[,c("species", unique(unlist(nn_vg[,1:2])))]  
nn_res <- NULL
#nn_res <- foreach(i = 7:nrow(nn_vg), .combine = rbind,
#.inorder = F, .errorhandling = "remove") %dopar%{
#if (!require("pacman")) install.packages("pacman")
#pacman::p_unload(pacman::p_loaded(), character.only = T)
#pacman::p_load(pkgs, character.only = T)
nn_ids <- get_vgr_id(run_vars, nn_vg) 
# nn_ids <- 1:nrow(nn_vg)
for(i in nn_ids){
  print(i)
  nn_res <- rbind(nn_res, 
                  pglsPhyloCor(pair = nn_vg[i, 1:2], 
                               data = nn_dat, tree = tree, 
                               log.vars = log.vars, 
                               pair.type = "nn", mgm_types = mgm_types, 
                               alco = F, overwrite = overwrite))
}

  #stopCluster(cl)  
nn_res <- nn_res[order(abs(nn_res$phylocor), decreasing = T),]
if(save){
  write.csv(nn_res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                          min.n, if(alco){"alco_"}, if(log){"_log"},
                          "_nn.csv", sep = ""),
            row.names = F)
}


# ---- pn-cc.list -------------------------------------------------------------------
# categorical - categorical edges (returns list)
cc_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(x == "c")}),]
cc_dat <- data[,c("species", unique(unlist(cc_vg[,1:2])))]  
cc_res.list <- vector("list", nrow(cc_vg))
cc_ids <- get_vgr_id(run_vars, cc_vg)
for(i in cc_ids){
  print(i)
  cc_res.list[[i]] <- pglsPhyloCor(pair = cc_vg[i, 1:2], data = cc_dat, 
                                   tree = tree, log.vars = NULL, pair.type = "cc",
                                   result = "row", overwrite = overwrite,
                                   mgm_types = mgm_types)
}
if(save){
  save(cc_res.list, file = paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                                 min.n, if(alco){"alco_"}, if(log){"_log"},"_cc.Rdata", 
                                 sep = ""))
}

# ---- pn-cc-load.list ----
load(file = paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn",
                  min.n, if(alco){"alco_"}, if(log){"_log"},"_cc.Rdata", sep = ""))

# ---- pn-cc ----
# categorical - categorical edges (list to data.frame)
cc_res <- NULL
for(i in 1:length(cc_res.list)){
  cc_res <- rbind(cc_res, get_cc.row(cc_res.list[[i]], cc_vg))
  print(i)
}
if(save){
  write.csv(cc_res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                          min.n, if(alco){"alco_"}, if(log){"_log"},"_cc.csv", sep = ""),
            row.names = F)}


# ---- pn-nc -------------------------------------------------------------------
# numeric - categorical edges
nc_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(c("c", "g") %in% x)}),]
nc_dat <- data[,c("species", unique(unlist(nc_vg[,1:2])))]
nc_res <- NULL
nc_ids <- get_vgr_id(run_vars, nc_vg) 
# nc_ids <- 1:nrow(nc_vg)
for(i in nc_ids){
  print(i)
  nc_res <- rbind(nc_res, 
                  pglsPhyloCor(pair = nc_vg[i, 1:2], data = nc_dat, 
                               tree = tree, log.vars = log.vars, 
                               pair.type = "nc", result = "row", 
                               mgm_types = mgm_types, alco = F, overwrite = overwrite))
}
if(save){write.csv(nc_res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                                 min.n, if(alco){"alco_"}, if(log){"_log"},"_nc.csv", sep = ""),
                   row.names = F)}


# ---- pn-read.type-res ----
nn_res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(alco){"alco_"}, if(log){"_log"},"_nn.csv", sep = ""))
nc_res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(alco){"alco_"}, if(log){"_log"},"_nc.csv", sep = ""))
cc_res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(alco){"alco_"}, if(log){"_log"},"_cc.csv", sep = ""))
# ---- pn-comp.res ----
res <- rbind(nn_res, nc_res, cc_res)
if(save){
  write.csv(res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                       min.n, if(alco){"alco_"}, if(log){"_log"},"_allDT.csv", sep = ""),
            row.names = F)
}


# ---- pn-read.type-res ----
res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                      min.n, if(alco){"alco_"}, if(log){"_log"},"_allDT.csv", sep = ""))

# ---- pn-run-netcarto ----
run_netcarto(res, edge_det = "p-value", alco, save = T)
run_netcarto(res, edge_det = "phylocor", alco, save = T)


