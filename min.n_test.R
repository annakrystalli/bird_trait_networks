rm(list = ls())
wkf = "goodmankruskal"
param = "min.n_test.R"
file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
source("~/Documents/workflows/bird_trait_networks/project_ui.R")


require(plotly)

min.n_df <- NULL

set.seed(1)
spl.n <- sort(sample(x = sort(unique(vg$n)), size = 200, replace = FALSE, 
                         prob = table(vg$n)/sum(table(vg$n))))

for(i in spl.n){
  vg.i <- vg[vg$n > i,]
  min.n_df <- rbind(min.n_df,
                    data.frame(min.n = i, pair.n = nrow(vg.i),
                               var.n = length(unique(unlist(vg.i[,1:2])))))
}



ay <- list(
  tickfont = list(color = "red"),
  overlaying = "y",
  side = "right"
)


plot_ly(min.n_df, x = min.n, y = pair.n, name = "pairwise n") %>%
  add_trace(x = min.n, y = var.n, name = "variable n", yaxis = "y2") %>%
  layout(title = "Effect of min.n on datapoint availability", yaxis2 = ay) 
    


i <- order(TDdf$Species)
plot_ly(TDdf, x = Species, y = Dplus, name = "mean TD", 
        mode = "markers", text = var) %>%
  add_trace(y = c(EDplus, EDplus), x= c(min(Species), max(Species)), mode = "lines") %>%
  add_trace(x = Species[i], y = EDplus[i] - 2 * sd.Dplus[i], name = "EDplus", 
            mode = "lines") %>%
  add_trace(x = Species[i], y = EDplus[i] + 2 * sd.Dplus[i], name = "EDplus", 
            mode = "lines") %>%
  layout(title = "Taxonomic diversity of variable subsamples according to data availability ")

    
    
  

  lines(x$Species[i], x$EDplus - 2 * x$sd.Dplus[i], ...)
  lines(x$Species[i], x$EDplus + 2 * x$sd.Dplus[i], ...)
