#' ---
#' output:
#' html_document:
#' self_contained: false
#' ---

## ---- load_dependence ----
lapply(c(pkgs, "qgraph"), require, character.only = TRUE, quietly = TRUE)

