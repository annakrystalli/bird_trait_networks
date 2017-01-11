## ---- source-params ----
if(exists("param")){
  source(paste(script.folder, "params/", param, sep = ""))
}
## ---- source-wkf ----
## load_rmacro
if(exists("wkf")){ 
  if(wkf == "rmacro"){
    eval(parse(text = getURL("https://raw.githubusercontent.com/annakrystalli/rmacroRDM/master/R/functions.R", ssl.verifypeer = FALSE)))
    eval(parse(text = getURL("https://raw.githubusercontent.com/annakrystalli/rmacroRDM/master/R/wideData_function.R", ssl.verifypeer = FALSE)))
  }
  ## load_phylocor
  if(wkf == "phylocor"){
    source(paste(script.folder, "PhyloCor functions.R", sep = ""))
  }
  ## load_phylonetworker
  if(wkf == "phylonetworker"){
    source(paste(script.folder, "PhyloNetworker_functions.R", sep = ""))
    source(paste(script.folder, "PhyloNetworker_setup.R", sep = ""))
  }
  # load_mcmc
  if(wkf == "mcmc"){
    source(paste(script.folder, "PhyloNetworker_functions.R", sep = ""))
    source(paste(script.folder, "PhyloNetworker_setup.R", sep = ""))
  }
  ## load_mgm
  if(wkf == "mgm"){
    source(paste(script.folder, "PhyloNetworker_functions.R", sep = ""))
    source(paste(script.folder, "PhyloNetworker_setup.R", sep = ""))
  }
  ## load_gk
  if(wkf == "goodmankruskal"){
    source(paste(script.folder, "PhyloNetworker_functions.R", sep = ""))
    source(paste(script.folder, "PhyloNetworker_setup.R", sep = ""))
  }
}
