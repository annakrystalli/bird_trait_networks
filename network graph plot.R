rm(list=ls())

source("~/Documents/workflows/bird_trait_networks/setup.R")


# PACKAGES & SETTINGS ###############################################################

library(igraph)

an.ID <- "100spp"
min.n <- 10
log <- T
cutoff <- 0.35
if(log){log.vars <- metadata$code[as.logical(metadata$log)]}else{log.vars <- ""}


# DATA ###############################################################################

load(file = paste(output.folder, "data/networks/", an.ID,"_net_mn", min.n, 
                  if(log){"_log"},".RData", sep = ""))

res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                            min.n, if(log){"_log"},".csv", sep = ""),
                 stringsAsFactors = F)


# process data ###############################################################################

edges <- res[abs(res$phylocor) > 0.35 & !is.na(res$phylocor),]

g.n <- as.vector(apply(edges[,c("var1", "var2")],1,FUN = unlist))
n.names <- unique(g.n)
g.code <- match(g.n, n.names)  

# node colours
role.levels <- levels(unique(net[[1]]$role))
role.codes <- 
  role.cols <- RColorBrewer::brewer.pal(7, "BrBG")
n.cols <- role.cols[match(as.character(net[[1]]$role)[match(n.names, net[[1]]$name)], role.levels)]

# edge colours
edge.cols <- rep("black", length(edges$phylocor))
edge.cols[edges$phylocor >= 0.1] <- "red"
edge.cols[edges$phylocor < -0.1] <- "blue"


# fit graph ###############################################################################

G <- graph(g.code, directed = FALSE)

# Assign attributes to the graph
G$name    <- "correlation network of numeric bird traits"

# Assign attributes to the graph's vertices
V(G)$name  <- 1:length(n.cols)
V(G)$color <- n.cols
E(G)$edge.color <- edge.cols


# Assign attributes to the edges
E(G)$weight <- edges$phylocor



# Plot the graph -- details in the "Drawing graphs" section of the igraph manual
png(filename = paste(output.folder, "figures/net.png", sep = ""),
    width =  1500, height = 1194, pointsize = 20)

par(mai = c(0,0.5,0.5,0))

plot(G, layout = layout.fruchterman.reingold, 
     main = G$name,
     vertex.label = V(G)$name,
     vertex.size = 5.5,
     vertex.color= V(G)$color,
     vertex.frame.color= "white",
     vertex.label.color = "black",
     vertex.label.family = "sans",
     edge.width=E(G)$weight,
     vertex.label.cex=0.65,
     edge.color= edge.cols)



legend("bottomleft", legend = role.levels[which(role.levels %in% as.character(net[[1]]$role))], 
       fill = role.cols[which(role.levels %in% as.character(net[[1]]$role))], 
       bty = "n", cex = 1)

dev.off()
