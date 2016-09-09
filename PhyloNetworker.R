# ---- init-an ----
source(paste(script.folder, "project_ui.R", sep = ""))

num.res <- NULL
for(i in 1:dim(num.vg)[1]){
  
  num.res <- rbind(num.res, pglsPhyloCor(x = num.vg[i, 1:2], data = num.dat, 
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

