#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---

# ---- install_dependencies ----
in.id <- sapply(pkgs, FUN = function(x) {require(x, character.only = TRUE)})
in.pkgs <- pkgs[!in.id]
  
#' install from CRAN
if(length(in.pkgs != 0)){
install.packages(in.pkgs)}

#' install from github
install_github('SachaEpskamp/qgraph')
