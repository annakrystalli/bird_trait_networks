rm(list = ls())

file_setup_path <- "file_setup.R"
source(file_setup_path)
require(rmarkdown)
require(knitr)
require(dplyr)

# ---- options ---- 
theme <- "journal"

# ---- render-site-ALL ----
files <- list.files("man/")[grep(".Rmd", list.files("man/"))]
files <- files[files != "Hierarchical Networks.Rmd"]
# files <- files[-1:-3]
for(file in files){
render(paste0("man/", file), output_format = "html_notebook",
       output_options = list(theme = theme),          
       output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))}

# ---- render-site-INDIVIDUAL ----
file <- "project_README.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))

  file <- "project_ui.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "file_system.Rmd"
  render(paste0("man/", file), output_format = "html_document", 
         output_options = list(theme = theme),
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".html", file)))
  
  file <- "PhyloNetworker_README.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "pglsPhyloCor_README.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "GoodmanKruskal.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "rmacroRDM.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "Hierarchical Networks.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "Network_Data_Availability_all.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "Network_Data_Availability_100spp.Rmd"
  render(paste0("man/", file), output_format = "html_notebook",
         output_options = list(theme = theme),          
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "pulling_repo.Rmd"
  render(paste0("man/", file), output_format = "html_document",
         output_options = list(theme = theme),
         output_file = paste0(script.folder,"docs/",gsub(".Rmd", ".html", file)))

  spin("mgm.R")
  file.copy(from = "mgm.html", to = "docs/mgm.html", overwrite = T)
  
  file <- "data_gap_eval.Rmd"
  render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
         output_options = list(theme = theme),  
         output_file =  paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  file <- "mgm_viz.Rmd"
  render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
         output_options = list(theme = theme),
         output_file =  paste0(script.folder,"docs/",gsub(".Rmd", ".nb.html", file)))
  
  # ---- render-results ----
  file <- "data_gap_eval.Rmd"
  render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
         output_dir = paste(output.folder, "Reports/Results/", sep = ""))
  
  file <- "Categorical_phylocors.Rmd"
  render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
         output_dir = paste(output.folder, "reports/results/", sep = ""))
  
  source(paste0(script.folder, "params/phylonetworker.R"))
  file <- "Network_viz.Rmd"
  render(paste("reports/results/", file, sep =""), output_format = "html_notebook",
         output_file = paste(output.folder,"Reports/Results/",gsub(".Rmd", paste0("_",an.ID), 
                                                                   file), ".nb.html", sep =""))

  