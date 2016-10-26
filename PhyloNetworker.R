rm(list = ls())
wkf = "phylonetworker"
param = "phylonetworker.R"
file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
source("~/Documents/workflows/bird_trait_networks/project_ui.R")
vg_dt <- get_vg_dt(vg, mgm_types)

# ---- nn ----
gg_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(x %in% c("g", "p"))}),]
gg_dat <- data[,c("species", unique(unlist(gg_vg[,1:2])))]  

num.res <- NULL
for(i in 1:dim(num.vg)[1]){
  
  num.res <- rbind(num.res, pglsPhyloCor(x = gg.vg[i, 1:2], data = gg.dat, 
                                         tree = tree, log.vars = log.vars, 
                                         datTypes = "nn"))
  print(i)
}


num.res <- num.res[order(abs(num.res$phylocor), decreasing = T),]

write.csv(num.res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_num.csv", sep = ""),
          row.names = F)


bin.res <- NULL
for(i in 1:dim(bin.vg)[1]){
  
  bin.res <- rbind(bin.res, pglsPhyloCor(x = bin.vg[i, 1:2], data = bin.dat, 
                                         tree = tree, log.vars = log.vars, 
                                         datTypes = "bb"))
  print(i)
}


bin.res <- bin.res[order(abs(bin.res$phylocor), decreasing = T),]

write.csv(bin.res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_bin.csv", sep = ""),
          row.names = F)



# ---- cc ----

cc_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(x == "c")}),]
cc_dat <- data[,c("species", unique(unlist(cc_vg[,1:2])))]  


cc.res <- pglsPhyloCor(pair = cc_vg[nrow(cc_vg)-20, 1:2], data = cc_dat, 
                       tree = tree, log.vars = NULL, datTypes = "cc", result = "row")




# ---- cat vs numeric (nc)----


vg <- orderby_vg_dt(vg, mgm_types)
vg_dt <- get_vg_dt(vg, mgm_types)
cg_vg <- vg[apply(vg_dt, 1, FUN = function(x) {all(c("c", "g") %in% x)}),]
cg_dat <- data[,c("species", unique(unlist(cg_vg[,1:2])))]

cg.res <- NULL

for(i in 7:dim(cg_vg)[1]){
  cg.res <- rbind(cg.res, pglsPhyloCor(pair = cg_vg[i, 1:2], data = cg_dat, 
                                         tree = tree, log.vars = log.vars, 
                                         datTypes = "nc", result = "row"))
  print(i)
}



write.csv(cg.res, paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                         min.n, if(log){"_log"},"_cg.csv", sep = ""),
          row.names = F)


dat <- cg.res[is.na(cg.res$error), ]

plot(cg.res$Dplus,  cg.res$lambda)
lm(lambda ~ Dplus, data = dat)
lm(phylocor ~ Dplus, data = dat)



### mixed models





#################################################################################
## NETWORK ANALYSIS ####################################################################
##################################################################################
res <- res[abs(res$phylocor) > cutoff & !is.na(res$phylocor),]

library(rnetcarto)
net.d <- res[!is.na(res$phylocor), c("var1", "var2", "phylocor")]
net.list <- as.list(net.d)
net <- netcarto(web = net.list, seed = 1)


write.csv(cbind(net[[1]], modularity = net[[2]]), 
          paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                if(log){"_log"},".csv", sep = ""),
          row.names = F)

save(net, file = paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                       if(log){"_log"},".RData", sep = ""))







#source('http://bioconductor.org/biocLite.R')
#biocLite ('RCytoscape')

