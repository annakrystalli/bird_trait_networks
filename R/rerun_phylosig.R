rm(list = ls())

# ---- pn-setup ----
if(exists("file_setup_path")){}else{
  file_setup_path <- "~/Documents/workflows/bird_trait_networks/file_setup.R"
  source(file_setup_path)}
wkf = "phylonetworker"
param = "phylonetworker.R"
source(paste0(script.folder,"project_ui.R"))

# ---- functions ----
res_amend.phylosig <- function(res, pair, data, tree, log.vars = NULL, pair.type, 
                               mgm_types = mgm_types, verbose = T){
  
  pair <- as.character(pair)
  t0 <- Sys.time()
  if(verbose){cat("processing pair:", pair, "; n:",  res[i, "n"], "\n")}
  if(!pair.type %in% c("nc", "nn")) {return(res)}
  
  # ---- set-up-vars ----
  data <- data[, c("species", gsub("log_", "", pair))]
  
  # ---- sub-complete-cases ----
  data <- data[complete.cases(data),]
  
  # ---- log-log.vars ----  
  for(var in pair){
    if(mgm_types[gsub("log_", "", var)] == "g"){
      if(validateLog(var = gsub("log_", "", var), log.vars, data)){
        data[, gsub("log_", "", var)] <- log(data[, gsub("log_", "", var)])
      }
    }
  }
  
  # ---- fit-amend ---- 
  amend <- try(get_physig(data, tree, pair, mgm_types = mgm_types, 
                          pair.type))
  if(class(amend) == "try-error"){
    amend <- NA
  }
  
  res[res$var1 == pair[1] & res$var2 == pair[2], 
      c("physig.var1", "physig.var2")] <- amend
  
  if(verbose){cat("time elapsed:", t0 - Sys.time(), "sec", "\n")
    cat("phylosig complete =", TRUE, "\n")}
  cat("***", "\n")
  
  # ---- return ----
  return(res)
  
}

get_physig <- function(data, tree, pair, mgm_types = mgm_types, 
                       pair.type) {
  
  # ---- configure-tree ----
  tree <- drop.tip(tree, setdiff(tree$tip.label, data$species))
  
  # ---- get-phylosig ----
  physig.var1 <- phylosig(tree, setNames(data[,gsub("log_", "", pair[1])], data$species),
                          method = "lambda")$lambda
  if(mgm_types[gsub("log_", "", pair[2])] == "g"){
    physig.var2 <- phylosig(tree, setNames(data[,gsub("log_", "", pair[2])], data$species),
                            method = "lambda")$lambda
  }else{
    physig.var2 <- NA  
  }
  out <- c(physig.var1, physig.var2)
  return(out)
}


# ---- pn-get-vg_dt ----
vg <- orderby_vg_dt(vg, mgm_types)
vg_dt <- get_vg_dt(vg, mgm_types)

res <- read.csv(paste(output.folder, "data/phylocors/", an.ID,"_phylocor_mn", 
                      min.n, if(log){"_log"},"_allDT.csv", sep = ""))
# res <- res[!apply(res, 1, function(x){any(pair[1] %in% x)}),]

i = 4
pair <- res[i, c("var1", "var2")]
pair.type <- res[i, "pair.type"]

k <- 1
for(i in which(res$pair.type %in% c("nn", "nc"))){
  cat("row.id: ", i, " (", k, " of ", sum(res$pair.type %in% c("nn", "nc")), ")", sep = "",
      "\n")
  res <- res_amend.phylosig(res, pair = res[i, c("var1", "var2")], data, tree, 
                            log.vars = log.vars, pair.type = res[i, "pair.type"], 
                          mgm_types = mgm_types, verbose = T)
  k <- k + 1
}

