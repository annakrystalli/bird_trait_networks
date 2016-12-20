
# ---- bb ----
# binary - binary edges

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
