rm(list = ls())

file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
source(file_setup_path)
require(rmarkdown)

file <- "mgm_viz.Rmd"
render(paste("reports/results/", file, sep =""),  output_format = "html_notebook",
       output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))

file <- "data_gap_eval.Rmd"
render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
output_dir = paste(output.folder, "Reports/Results/", sep = ""))

file <- "Categorical_phylocors.Rmd"
render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
       output_dir = paste(output.folder, "Reports/Results/", sep = ""))

source(paste0(script.folder, "params/phylonetworker.R"))
file <- "Network_viz.Rmd"
render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
       output_file = paste(output.folder,"reports/results/",gsub(".Rmd", paste0("_",an.ID), 
                                                         file), ".nb.html", sep =""))


# render site
# 
files <- c("PhyloNetworker_README.Rmd", "pglsPhyloCor_README.Rmd")
for(file in files){
render(paste0("man/", file), output_format = "html_notebook",
       output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))}

file <- "project_README.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))

  file <- "project_ui.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "GoodmanKruskal.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "rmacroRDM.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  

  spin("mgm.R")
  file.copy(from = "mgm.html", to = "docs/mgm.html", overwrite = T)
  
  file <- "GoodmanKruskal.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  