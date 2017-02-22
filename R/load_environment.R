## ---- source-params ----
if(exists("param")){
  source(paste(script.folder, "params/", param, sep = ""))
}
## ---- source-wkf ----
## load_rmacro
if(exists("wkf") & exists("param")){ 
  if(wkf == "rmacro"){
    eval(
      parse(text = 
              getURL("https://raw.githubusercontent.com/annakrystalli/rmacroRDM/master/R/functions.R", 
                     ssl.verifypeer = FALSE)
      )
    )
    sr_configurator <- parse(text = 
                               getURL("https://raw.githubusercontent.com/annakrystalli/rmacroRDM/master/R/sys.ref_configurator.R", 
                                      ssl.verifypeer = FALSE))
  }
  
  ## load_phylonetworker
  if(wkf == "phylonetworker"){
    source(paste(script.folder, "R/PhyloNetworker_functions.R", sep = ""))
    source(paste(script.folder, "R/PhyloNetworker_setup.R", sep = ""))
  }
  ## load_mgm
  if(wkf == "mgm"){
    source(paste(script.folder, "R/PhyloNetworker_functions.R", sep = ""))
    source(paste(script.folder, "R/PhyloNetworker_setup.R", sep = ""))
  }
  ## load_gk
  if(wkf == "goodmankruskal"){
    source(paste(script.folder, "R/PhyloNetworker_functions.R", sep = ""))
    source(paste(script.folder, "R/PhyloNetworker_setup.R", sep = ""))
  }
}
