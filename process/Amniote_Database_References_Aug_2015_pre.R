## ---- prR1-subset ----
data <- data[data$class == "Aves",]

## ---- prR1-derive_repro.age.diff ----
data$repro.age.diff <- paste(data$female_maturity_d, data$male_maturity_d, sep = "; ")
data$repro.age.diff[is.na(data$female_maturity_d) | is.na(data$male_maturity_d)] <- NA
data <- data[, !names(data) %in% var.omit]

add_keep.vars <- "repro.age.diff"
