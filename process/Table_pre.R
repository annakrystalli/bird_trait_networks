var.col <- sr$vnames[sr$vnames$code == "var", dcode]
data[data[, var.col] == "Habitat Foraging", var.col] <- "Habitat foraging"



